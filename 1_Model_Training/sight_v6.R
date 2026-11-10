# === SIGHT v6.5 (Fast Screening via 1-Fold, Full 10-Fold on Best Model) ===
# 1. Screen all models on 1 validation holdout across thresholds 0.05-0.95.
# 2. Pick the best model/params/threshold optimizing F2 (favoring sensitivity).
# 3. Run full 10-fold CV only on the single winning configuration.

library(tidyverse)
library(tidymodels)
library(workflowsets)
library(xgboost)
library(lightgbm)
library(ranger)
library(glmnet)
library(vip)
library(doParallel)
library(ggplot2)
library(yardstick)
library(patchwork)
library(bonsai)

# =======================================================
# === 1. LOAD AND PREPARE DATA ===
# =======================================================
tryCatch({
  df <- read_csv("combined_cleaned_data_for_training.csv")
}, error = function(e) {
  stop("Failed to load data file.")
})

POPULATION_COLUMN_NAME <- "source"
if (!POPULATION_COLUMN_NAME %in% names(df)) {
  df <- df %>% mutate(!!POPULATION_COLUMN_NAME := "Unknown")
}

df_prepared <- df %>%
  mutate(
    scd_status = case_when(
      genotype %in% c("SS", "SC", "Sbeta0", "Sbeta+", "SS- Lepore", "SS/HPFH", "SF", "S-BETA") ~ "SCD",
      genotype %in% c("AA", "Trait", "AS", "CC", "AC") ~ "No_SCD",
      TRUE ~ NA_character_
    ),
    scd_status = factor(scd_status, levels = c("SCD", "No_SCD")) 
  ) %>%
  select(patient_mrn, genotype, scd_status, !!sym(POPULATION_COLUMN_NAME),
         hematocrit, mchc, mcv, rbc_count, total_hgb, rdw) %>%
  drop_na()

# =======================================================
# === 2. SPLITTING & FAST SCREENING RESAMPLE ===
# =======================================================
set.seed(123)
data_split <- initial_split(df_prepared, prop = 0.80, strata = scd_status)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Scale weights (Control / Case)
case_count <- sum(train_data$scd_status == "SCD")
control_count <- sum(train_data$scd_status == "No_SCD")
scale_pos_weight <- control_count / case_count

# 1-Fold Resample: Carve out a single validation split (80/20) from train_data
set.seed(42)
screen_resample <- validation_split(train_data, prop = 0.80, strata = scd_status)

model_recipe <- recipe(scd_status ~ ., data = train_data) %>%
  update_role(patient_mrn, genotype, !!sym(POPULATION_COLUMN_NAME), new_role = "ID") %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

# =======================================================
# === 3. MODEL SPECS ===
# =======================================================
xgb_spec <- boost_tree(
  trees = tune(), tree_depth = tune(), min_n = tune(), loss_reduction = tune(), 
  sample_size = tune(), learn_rate = tune(), mtry = tune()
) %>% set_engine("xgboost", scale_pos_weight = scale_pos_weight, counts = FALSE) %>% set_mode("classification")

lgbm_spec <- boost_tree(
  trees = tune(), tree_depth = tune(), min_n = tune(), loss_reduction = tune(), 
  sample_size = tune(), learn_rate = tune(), mtry = tune()
) %>% set_engine("lightgbm", objective = "binary", scale_pos_weight = scale_pos_weight, counts = FALSE) %>% set_mode("classification")

rf_spec <- rand_forest(trees = 500, mtry = tune(), min_n = tune()) %>% 
  set_engine("ranger", importance = "permutation") %>% set_mode("classification")

lr_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>% 
  set_engine("glmnet") %>% set_mode("classification")

# =======================================================
# === 4. WORKFLOW & FAST SCREENING (1 FOLD ONLY) ===
# =======================================================
model_wfs <- workflow_set(
  pre = list(recipe = model_recipe),
  models = list(xgboost = xgb_spec, lightgbm = lgbm_spec, rf = rf_spec, logreg = lr_spec),
  cross = TRUE
)

TUNE_GRID_SIZE <- 20
boost_grid <- grid_space_filling(
  trees(), tree_depth(), min_n(), loss_reduction(), sample_size = sample_prop(), 
  learn_rate(), mtry = mtry_prop(range = c(0.5, 1.0)), size = TUNE_GRID_SIZE
)
rf_params <- parameters(min_n(), mtry()) %>% finalize(train_data)
rf_grid <- grid_space_filling(rf_params, size = TUNE_GRID_SIZE)
lr_grid <- grid_space_filling(penalty(), mixture(), size = TUNE_GRID_SIZE)

model_wfs <- model_wfs %>%
  option_add(grid = boost_grid, id = "recipe_xgboost") %>%
  option_add(grid = boost_grid, id = "recipe_lightgbm") %>%
  option_add(grid = rf_grid, id = "recipe_rf") %>%
  option_add(grid = lr_grid, id = "recipe_logreg")

# Parallel Execution for Fast 1-Fold Search
cores <- parallel::detectCores(logical = FALSE)
cl <- makePSOCKcluster(max(1, cores - 1))
clusterExport(cl, "scale_pos_weight")
registerDoParallel(cl)

cat("--- Fast Screening on Single Validation Holdout (1 Fold)... ---\n")
screening_results <- model_wfs %>%
  workflow_map("tune_grid", resamples = screen_resample, 
               metrics = metric_set(roc_auc),
               control = control_grid(save_pred = TRUE, verbose = FALSE))

# =======================================================
# === 5. SCREEN THRESHOLDS (0.05 to 0.95) & PICK WINNER ===
# =======================================================
cat("--- Searching Thresholds (0.05 - 0.95) Favoring Sensitivity (F2)... ---\n")

threshold_candidates <- seq(0.05, 0.95, by = 0.10)

evaluate_thresholds_for_tune <- function(wf_set_res, wf_id, thresholds) {
  preds <- collect_predictions(extract_workflow_set_result(wf_set_res, wf_id))
  
  # Group by each tuned hyperparameter config (.config)
  preds %>%
    group_by(.config) %>%
    group_modify(~ {
      df_grp <- .x
      map_dfr(thresholds, function(th) {
        df_eval <- df_grp %>%
          mutate(pred_class = factor(ifelse(.pred_SCD >= th, "SCD", "No_SCD"), levels = c("SCD", "No_SCD")))
        
        tibble(
          threshold   = th,
          sensitivity = sens(df_eval, truth = scd_status, estimate = pred_class)$.estimate,
          specificity = spec(df_eval, truth = scd_status, estimate = pred_class)$.estimate,
          f2_score    = f_meas(df_eval, truth = scd_status, estimate = pred_class, beta = 2)$.estimate,
          accuracy    = accuracy(df_eval, truth = scd_status, estimate = pred_class)$.estimate
        )
      })
    }) %>%
    ungroup() %>%
    mutate(wflow_id = wf_id)
}

# Run threshold evaluation over all screened models
all_threshold_perf <- map_dfr(screening_results$wflow_id, ~ {
  evaluate_thresholds_for_tune(screening_results, .x, threshold_candidates)
})

# Pick the single combination that maximizes F2 score, breaking ties by sensitivity
best_candidate <- all_threshold_perf %>%
  arrange(desc(f2_score), desc(sensitivity)) %>%
  slice(1)

BEST_WFLOW_ID   <- best_candidate$wflow_id
BEST_CONFIG     <- best_candidate$.config
FINAL_THRESHOLD <- best_candidate$threshold

cat(sprintf("\n🏆 Top Screening Configuration:\n Model: %s\n Threshold: %.2f\n Validation F2: %.3f\n Validation Sens: %.3f\n Validation Spec: %.3f\n",
            BEST_WFLOW_ID, FINAL_THRESHOLD, best_candidate$f2_score, best_candidate$sensitivity, best_candidate$specificity))

# Extract the exact hyperparameter values corresponding to the best configuration
wf_result <- extract_workflow_set_result(screening_results, BEST_WFLOW_ID)
best_params <- wf_result %>%
  collect_metrics() %>%
  filter(.config == BEST_CONFIG) %>%
  distinct(across(-c(.metric, .estimator, mean, n, std_err, .config)))

finalized_wf <- finalize_workflow(extract_workflow(screening_results, BEST_WFLOW_ID), best_params)

# =======================================================
# === 6. RUN FULL 10-FOLD CV ONLY ON THE WINNING MODEL ===
# =======================================================
cat("\n--- Running Full 10-Fold CV ONLY on the Winning Model... ---\n")
set.seed(42)
cv_folds <- vfold_cv(train_data, v = 10, strata = scd_status)

cv_results_winning <- fit_resamples(
  finalized_wf,
  resamples = cv_folds,
  metrics = metric_set(roc_auc, accuracy, sens, spec, ppv),
  control = control_resamples(save_pred = TRUE)
)

stopCluster(cl)

# Display 10-fold CV results for the selected threshold
cv_preds <- collect_predictions(cv_results_winning) %>%
  mutate(pred_class = factor(ifelse(.pred_SCD >= FINAL_THRESHOLD, "SCD", "No_SCD"), levels = c("SCD", "No_SCD")))

cv_threshold_metrics <- metric_set(sens, spec, ppv, npv, accuracy)(cv_preds, truth = scd_status, estimate = pred_class)
cat(sprintf("\n10-Fold CV Performance of Winning Model at Threshold %.2f:\n", FINAL_THRESHOLD))
print(cv_threshold_metrics)

# =======================================================
# === 7. FINAL FIT ON FULL TRAIN, EVALUATE ON TEST ===
# =======================================================
cat("\n--- Fitting on Entire Training Set and Testing... ---\n")
set.seed(123)
final_fit <- last_fit(finalized_wf, data_split)
test_predictions <- collect_predictions(final_fit)

final_results <- test_predictions %>%
  bind_cols(test_data %>% select(!!sym(POPULATION_COLUMN_NAME))) %>%
  mutate(.final_pred_class = factor(ifelse(.pred_SCD >= FINAL_THRESHOLD, "SCD", "No_SCD"), levels = c("SCD", "No_SCD")))

cat("\nFinal Test Metrics:\n")
test_perf <- metric_set(sens, spec, ppv, npv, accuracy)(final_results, truth = scd_status, estimate = .final_pred_class)
print(test_perf)

# =======================================================
# === 8. EXPORT FIGURES ===
# =======================================================
dir.create("publication_figures", showWarnings = FALSE)
theme_set(theme_bw(base_size = 14))

# -- Fig 1: Confusion Matrix --
p1 <- conf_mat(final_results, truth = scd_status, estimate = .final_pred_class) %>% autoplot(type = "heatmap") +
  scale_fill_gradient(low = "#F0F0F0", high = "#3B82F6") +
  labs(title = "Confusion Matrix", subtitle = paste("Optimized Threshold:", FINAL_THRESHOLD))
ggsave("publication_figures/fig1_confusion_matrix.png", p1, width = 6, height = 5)

# -- Fig 3: Subgroup Metrics --
subgroup_mets <- final_results %>% group_by(!!sym(POPULATION_COLUMN_NAME)) %>% 
  metric_set(accuracy, sens, spec)(truth = scd_status, estimate = .final_pred_class)

p3 <- ggplot(subgroup_mets, aes(x = !!sym(POPULATION_COLUMN_NAME), y = .estimate, fill = .metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = scales::percent(.estimate, 0.1)), position = position_dodge(0.9), vjust = -0.2, size = 3) +
  scale_fill_brewer(palette = "Set2") + labs(title = "Performance by Population")
ggsave("publication_figures/fig3_subgroup_metrics.png", p3, width = 9, height = 6)

# -- Fig 5: ROC Curve --
p5 <- roc_curve(final_results, truth = scd_status, .pred_SCD) %>% autoplot() +
  labs(title = "ROC Curve", subtitle = paste("AUROC:", round(roc_auc(final_results, scd_status, .pred_SCD)$.estimate, 3)))
ggsave("publication_figures/fig5_roc_curve.png", p5, width = 6, height = 6)

# -- Fig 7: Probability Distribution --
p7 <- ggplot(final_results, aes(x = .pred_SCD, fill = scd_status)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = FINAL_THRESHOLD, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = FINAL_THRESHOLD + 0.03, y = 0, label = paste("Cutoff:", FINAL_THRESHOLD), 
           color = "red", angle = 90, vjust = -0.5, hjust = 0) +
  scale_fill_manual(values = c("SCD" = "red", "No_SCD" = "blue")) +
  labs(title = "Predicted Probability Distribution",
       subtitle = paste("Optimized Threshold:", FINAL_THRESHOLD),
       x = "Predicted Probability (SCD)", y = "Density") +
  theme_minimal()
ggsave("publication_figures/fig7_prob_density.png", p7, width = 7, height = 5)

# ==============================================================================
# === 9. LANCET HAEMATOLOGY PUBLICATION SVG EXPORT MODULE (SIGHT) ===
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(scales)
library(patchwork)
library(svglite)
library(shapviz)

sight_svg_dir <- "publication_figures_SIGHT_svg"
dir.create(sight_svg_dir, showWarnings = FALSE, recursive = TRUE)

# Official Lancet Color Palette (ggsci 'lanonc')
lancet_pal <- c(
  "Navy"    = "#00468B",
  "Crimson" = "#ED0000",
  "Teal"    = "#42B540",
  "Cyan"    = "#0099B4",
  "Purple"  = "#925E9F",
  "Peach"   = "#FDAF91",
  "Red"     = "#AD002A",
  "Slate"   = "#ADB6B6",
  "Charcoal"= "#1B1919"
)

# Unified Lancet Theme with Times New Roman
theme_lancet_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "Times New Roman") +
    theme(
      plot.title       = element_text(family = "Times New Roman", size = base_size + 2, face = "bold", color = "#1B1919", hjust = 0),
      plot.subtitle    = element_text(family = "Times New Roman", size = base_size, color = "#424242", hjust = 0, margin = margin(b = 8)),
      axis.title       = element_text(family = "Times New Roman", size = base_size, face = "bold", color = "#1B1919"),
      axis.text        = element_text(family = "Times New Roman", size = base_size - 1, color = "#1B1919"),
      legend.title     = element_text(family = "Times New Roman", size = base_size - 0.5, face = "bold", color = "#1B1919"),
      legend.text      = element_text(family = "Times New Roman", size = base_size - 1, color = "#1B1919"),
      strip.text       = element_text(family = "Times New Roman", size = base_size, face = "bold", color = "#1B1919"),
      strip.background = element_rect(fill = "#F4F6F8", color = "#CCCCCC"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EAEAEA", linewidth = 0.35),
      panel.border     = element_rect(color = "#1B1919", fill = NA, linewidth = 0.6),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.tag         = element_text(family = "Times New Roman", size = base_size + 4, face = "bold", color = "#1B1919")
    )
}

# Standard English Feature Name Dictionary
clean_sight_features <- function(x) {
  recode(x,
         "hematocrit" = "Hematocrit (%)",
         "rdw"        = "RDW (%)",
         "total_hgb"  = "Total Hemoglobin (g/dL)",
         "mcv"        = "MCV (fL)",
         "mchc"       = "MCHC (g/dL)",
         "rbc_count"  = "RBC Count (10^6/uL)",
         .default     = x
  )
}

# ------------------------------------------------------------------------------
# Panel 2A: SIGHT Confusion Matrix (Custom ggplot tile for full label control)
# ------------------------------------------------------------------------------
cm_sight_df <- final_results %>%
  mutate(
    Truth = factor(ifelse(scd_status == "SCD", "SCD", "No SCD"), levels = c("SCD", "No SCD")),
    Prediction = factor(ifelse(.final_pred_class == "SCD", "SCD", "No SCD"), levels = c("No SCD", "SCD"))
  ) %>%
  count(Truth, Prediction, .drop = FALSE) %>%
  group_by(Truth) %>%
  mutate(
    pct = n / sum(n) * 100,
    cell_label = sprintf("%s\n(%.1f%%)", comma(n), pct),
    text_col = ifelse(n > max(n) * 0.45, "white", "#1B1919")
  ) %>%
  ungroup()

fig2a_cm <- ggplot(cm_sight_df, aes(x = Truth, y = Prediction, fill = n)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_label, color = text_col),
            family = "Times New Roman", fontface = "bold", size = 4) +
  scale_color_identity() +
  scale_fill_gradient(low = "#EBF2F7", high = "#00468B", name = "Patient Count") +
  labs(
    title = "Confusion Matrix (Validation Cohort)",
    subtitle = sprintf("Optimized Classification Threshold: %.2f", FINAL_THRESHOLD),
    x = "Reference Standard (True Status)",
    y = "Predicted Status"
  ) +
  theme_lancet_pub() +
  theme(panel.grid = element_blank(), legend.position = "none")

ggsave(file.path(sight_svg_dir, "Fig2A_SIGHT_Confusion_Matrix.svg"), fig2a_cm, width = 5.5, height = 4.8, device = svglite)

# ------------------------------------------------------------------------------
# Panel 2B: SIGHT ROC Curve
# ------------------------------------------------------------------------------
roc_sight_df <- roc_curve(final_results, truth = scd_status, .pred_SCD)
auroc_sight_val <- roc_auc(final_results, truth = scd_status, .pred_SCD)$.estimate

fig2b_roc <- ggplot(roc_sight_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_abline(lty = 2, color = "#ADB6B6", linewidth = 0.7) +
  geom_path(color = "#00468B", linewidth = 1.2) +
  annotate("label", x = 0.62, y = 0.20,
           label = sprintf("AUROC = %.4f", auroc_sight_val),
           family = "Times New Roman", fontface = "bold", size = 3.8,
           fill = "#F4F6F8", color = "#00468B", label.size = 0.4) +
  scale_x_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  scale_y_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  labs(
    title = "Receiver Operating Characteristic",
    subtitle = "SIGHT Frontline CBC Screening Model",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) +
  theme_lancet_pub()

ggsave(file.path(sight_svg_dir, "Fig2B_SIGHT_ROC_Curve.svg"), fig2b_roc, width = 5.5, height = 4.8, device = svglite)

# ------------------------------------------------------------------------------
# Panel 2C: Performance by Population Subgroup
# ------------------------------------------------------------------------------
subgroup_clean_df <- final_results %>%
  group_by(!!sym(POPULATION_COLUMN_NAME)) %>%
  metric_set(accuracy, sens, spec)(truth = scd_status, estimate = .final_pred_class) %>%
  ungroup() %>%
  filter(!is.na(.estimate)) %>%
  mutate(
    Metric = factor(
      recode(.metric,
             "accuracy" = "Accuracy",
             "sens"     = "Sensitivity",
             "spec"     = "Specificity"),
      levels = c("Accuracy", "Sensitivity", "Specificity")
    ),
    Cohort = as.character(!!sym(POPULATION_COLUMN_NAME))
  )

fig2c_subgroups <- ggplot(subgroup_clean_df, aes(x = Cohort, y = .estimate, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.72, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", .estimate * 100)),
            position = position_dodge(width = 0.82), vjust = -0.45,
            family = "Times New Roman", fontface = "bold", size = 3.0, color = "#1B1919") +
  scale_fill_manual(
    values = c("Accuracy" = "#00468B", "Sensitivity" = "#ED0000", "Specificity" = "#0099B4"),
    name = "Performance Metric"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.10), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Diagnostic Performance by Cohort",
    subtitle = "Cross-Population Generalizability Evaluation",
    x = "Clinical / Geographic Cohort",
    y = "Performance Estimate (%)"
  ) +
  theme_lancet_pub() +
  theme(legend.position = "top")

ggsave(file.path(sight_svg_dir, "Fig2C_SIGHT_Subgroup_Performance.svg"), fig2c_subgroups, width = 7.0, height = 4.8, device = svglite)

# ------------------------------------------------------------------------------
# Panel 2D: SIGHT SHAP Feature Importance
# ------------------------------------------------------------------------------
fitted_sight_wf <- extract_workflow(final_fit)
saveRDS(fitted_sight_wf, "scd_final_lightgbm_model.rds")
saveRDS(fitted_sight_wf, "scd_final_model.rds")
sight_engine_fit <- extract_fit_parsnip(fitted_sight_wf)$fit
sight_baked_x <- extract_recipe(fitted_sight_wf) %>%
  bake(new_data = test_data, all_predictors()) %>%
  as.matrix()

sight_shap_df <- tryCatch({
  sv_sight <- shapviz(sight_engine_fit, X_pred = sight_baked_x, X = as.data.frame(sight_baked_x))
  if (is.list(sv_sight) && !inherits(sv_sight, "shapviz")) sv_sight <- sv_sight[[1]]
  imp_mat <- abs(get_shap_values(sv_sight))
  tibble(
    Feature = clean_sight_features(colnames(imp_mat)),
    Importance = colMeans(imp_mat, na.rm = TRUE)
  ) %>%
    arrange(desc(Importance))
}, error = function(e) {
  # Fallback to model-based importance if tree SHAP structure varies
  vi_df <- vip::vi(extract_fit_parsnip(fitted_sight_wf))
  tibble(
    Feature = clean_sight_features(vi_df$Variable),
    Importance = vi_df$Importance
  ) %>%
    arrange(desc(Importance))
})

fig2d_shap <- ggplot(sight_shap_df, aes(x = Importance, y = fct_reorder(Feature, Importance))) +
  geom_col(fill = "#00468B", width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", Importance)),
            hjust = -0.15, family = "Times New Roman", fontface = "bold", size = 3.2, color = "#1B1919") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "SHAP Predictor Importance",
    subtitle = "Mean Absolute SHAP Contribution to SCD Status",
    x = "Mean |SHAP Value| (Impact on Model Output)",
    y = "CBC Parameter"
  ) +
  theme_lancet_pub()

ggsave(file.path(sight_svg_dir, "Fig2D_SIGHT_SHAP_Importance.svg"), fig2d_shap, width = 6.0, height = 4.8, device = svglite)

# ------------------------------------------------------------------------------
# Combine Panels into Main Text Figure 2 (2x2 Multi-Panel Layout)
# ------------------------------------------------------------------------------
fig2_combined <- (fig2a_cm | fig2b_roc) / (fig2c_subgroups | fig2d_shap) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 15))

ggsave(file.path(sight_svg_dir, "Figure2_SIGHT_Main_4Panel.svg"), fig2_combined, width = 13.0, height = 10.5, device = svglite)

# ------------------------------------------------------------------------------
# Supplementary Figure S2: Cohort CBC Parameter Distributions
# ------------------------------------------------------------------------------
cbc_long_df <- df_prepared %>%
  mutate(
    Status = factor(ifelse(scd_status == "SCD", "SCD", "No SCD"), levels = c("SCD", "No SCD")),
    Cohort = as.character(!!sym(POPULATION_COLUMN_NAME))
  ) %>%
  pivot_longer(
    cols = c(hematocrit, mchc, mcv, rbc_count, rdw, total_hgb),
    names_to = "Parameter",
    values_to = "Value"
  ) %>%
  mutate(Parameter = clean_sight_features(Parameter))

supp_fig_s2 <- ggplot(cbc_long_df, aes(x = Cohort, y = Value, fill = Status)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.35, alpha = 0.88, width = 0.65, color = "#1B1919", linewidth = 0.35) +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("SCD" = "#ED0000", "No SCD" = "#00468B"), name = "Disease Status") +
  labs(
    title = "Distribution of Complete Blood Count Parameters Across Geographic Cohorts",
    subtitle = "Comparison of Baseline Hematological Biomarkers by Clinical Source and SCD Status",
    x = "Cohort Source",
    y = "Measured Parameter Value"
  ) +
  theme_lancet_pub() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(file.path(sight_svg_dir, "Supp_FigS2_SIGHT_Cohort_CBC_Distributions.svg"), supp_fig_s2, width = 12.0, height = 8.0, device = svglite)

# ------------------------------------------------------------------------------
# Supplementary Figure S3: Probability Density & Threshold Optimization Curve
# ------------------------------------------------------------------------------
prob_dens_df <- final_results %>%
  mutate(Status = factor(ifelse(scd_status == "SCD", "SCD", "No SCD"), levels = c("SCD", "No SCD")))

supp_s3a_density <- ggplot(prob_dens_df, aes(x = .pred_SCD, fill = Status, color = Status)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_vline(xintercept = FINAL_THRESHOLD, linetype = "dashed", color = "#AD002A", linewidth = 0.9) +
  annotate("label", x = FINAL_THRESHOLD + 0.14, y = 15,
           label = sprintf("Decision Threshold: %.2f", FINAL_THRESHOLD),
           family = "Times New Roman", fontface = "bold", size = 3.4,
           fill = "white", color = "#AD002A", label.size = 0.3) +
  scale_fill_manual(values = c("SCD" = "#ED0000", "No SCD" = "#00468B"), name = "True Status") +
  scale_color_manual(values = c("SCD" = "#ED0000", "No SCD" = "#00468B"), name = "True Status") +
  scale_x_continuous(labels = number_format(accuracy = 0.01), limits = c(0, 1)) +
  labs(
    title = "Predicted Probability Distribution",
    subtitle = "Bimodal Separation of SCD Cases and Controls on Validation Set",
    x = "Predicted Probability of SCD",
    y = "Kernel Density"
  ) +
  theme_lancet_pub() +
  theme(legend.position = "top")

thresh_curve_df <- all_threshold_perf %>%
  filter(wflow_id == BEST_WFLOW_ID, .config == BEST_CONFIG) %>%
  select(threshold, Sensitivity = sensitivity, Specificity = specificity, `F2-Score` = f2_score, Accuracy = accuracy) %>%
  pivot_longer(cols = -threshold, names_to = "Metric", values_to = "Score") %>%
  mutate(Metric = factor(Metric, levels = c("F2-Score", "Sensitivity", "Specificity", "Accuracy")))

supp_s3b_thresh <- ggplot(thresh_curve_df, aes(x = threshold, y = Score, color = Metric)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  geom_vline(xintercept = FINAL_THRESHOLD, linetype = "dashed", color = "#AD002A", linewidth = 0.8) +
  scale_color_manual(
    values = c("F2-Score" = "#ED0000", "Sensitivity" = "#00468B", "Specificity" = "#0099B4", "Accuracy" = "#42B540"),
    name = "Metric"
  ) +
  scale_x_continuous(breaks = seq(0.05, 0.95, by = 0.10), labels = number_format(accuracy = 0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.5, 1.0)) +
  labs(
    title = "Classification Threshold Optimization",
    subtitle = sprintf("Maximizing Sensitivity-Weighted F2-Score (Optimal Cutoff = %.2f)", FINAL_THRESHOLD),
    x = "Candidate Probability Threshold",
    y = "Metric Value (%)"
  ) +
  theme_lancet_pub() +
  theme(legend.position = "top")

supp_fig_s3 <- (supp_s3a_density | supp_s3b_thresh) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 14))

ggsave(file.path(sight_svg_dir, "Supp_FigS3_SIGHT_Prob_Density_and_Thresholds.svg"), supp_fig_s3, width = 12.5, height = 5.5, device = svglite)
ggsave(file.path(sight_svg_dir, "Supp_FigS3A_SIGHT_Prob_Density.svg"), supp_s3a_density, width = 6.2, height = 5.0, device = svglite)
ggsave(file.path(sight_svg_dir, "Supp_FigS3B_SIGHT_Threshold_Optimization.svg"), supp_s3b_thresh, width = 6.2, height = 5.0, device = svglite)

cat(sprintf("\n✅ All SIGHT Lancet SVGs exported successfully to '%s/'\n", sight_svg_dir))
# ==============================================================================
# === 10. LANCET HAEMATOLOGY & TRIPOD+AI TABLE EXPORT MODULE (SIGHT) ===
# ==============================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(pROC)

sight_table_dir <- "publication_tables_SIGHT"
dir.create(sight_table_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: Format Median (IQR) cleanly for Lancet style
fmt_med_iqr <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) == 0) return("NA")
  q <- quantile(x_clean, probs = c(0.25, 0.50, 0.75))
  sprintf("%.1f (%.1f–%.1f)", q[2], q[1], q[3])
}

# Helper: Clopper-Pearson Exact 95% CI for Proportions (Acc, Sens, Spec, PPV, NPV)
calc_binom_ci <- function(k, n) {
  if (is.na(k) || is.na(n) || n == 0) return("NA")
  bt <- binom.test(k, n)
  sprintf("%.1f%% (%.1f–%.1f)", (k / n) * 100, bt$conf.int[1] * 100, bt$conf.int[2] * 100)
}

# Helper: DeLong 95% CI for AUROC
calc_auroc_ci <- function(truth_vec, prob_vec) {
  if (length(unique(truth_vec[!is.na(truth_vec)])) < 2) return("NA (Single Class)")
  tryCatch({
    roc_obj <- pROC::roc(response = truth_vec, predictor = prob_vec, levels = c("No_SCD", "SCD"), direction = "<", quiet = TRUE)
    ci_vals <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    sprintf("%.4f (%.4f–%.4f)", ci_vals[2], ci_vals[1], ci_vals[3])
  }, error = function(e) "NA")
}

# ------------------------------------------------------------------------------
# 1. Table 1 (Part A): SIGHT Baseline Characteristics & Missing Data Audit
# ------------------------------------------------------------------------------
# Audit missing data from raw 'df' before drop_na() (Required by TRIPOD+AI)
missing_audit <- df %>%
  summarise(
    Cohort = "Missing Data Count (Pre-Filtering, N)",
    `Total N` = as.character(nrow(df) - nrow(df_prepared)),
    `SCD Cases, n (%)` = "—",
    `No SCD Controls, n (%)` = "—",
    `Hematocrit (%)` = as.character(sum(is.na(hematocrit))),
    `Total Hemoglobin (g/dL)` = as.character(sum(is.na(total_hgb))),
    `MCV (fL)` = as.character(sum(is.na(mcv))),
    `MCHC (g/dL)` = as.character(sum(is.na(mchc))),
    `RDW (%)` = as.character(sum(is.na(rdw))),
    `RBC Count (10^6/uL)` = as.character(sum(is.na(rbc_count)))
  )

summarise_sight_cohort <- function(data_sub, cohort_label) {
  n_tot <- nrow(data_sub)
  n_scd <- sum(data_sub$scd_status == "SCD", na.rm = TRUE)
  n_ctl <- sum(data_sub$scd_status == "No_SCD", na.rm = TRUE)
  
  tibble(
    Cohort = cohort_label,
    `Total N` = format(n_tot, big.mark = ","),
    `SCD Cases, n (%)` = sprintf("%s (%.1f%%)", format(n_scd, big.mark = ","), (n_scd / n_tot) * 100),
    `No SCD Controls, n (%)` = sprintf("%s (%.1f%%)", format(n_ctl, big.mark = ","), (n_ctl / n_tot) * 100),
    `Hematocrit (%)` = fmt_med_iqr(data_sub$hematocrit),
    `Total Hemoglobin (g/dL)` = fmt_med_iqr(data_sub$total_hgb),
    `MCV (fL)` = fmt_med_iqr(data_sub$mcv),
    `MCHC (g/dL)` = fmt_med_iqr(data_sub$mchc),
    `RDW (%)` = fmt_med_iqr(data_sub$rdw),
    `RBC Count (10^6/uL)` = fmt_med_iqr(data_sub$rbc_count)
  )
}

table1_sight_rows <- bind_rows(
  summarise_sight_cohort(df_prepared, "All SIGHT Participants (Combined)"),
  summarise_sight_cohort(train_data, "Training Set (80% Split)"),
  summarise_sight_cohort(test_data, "Validation Test Set (20% Split)"),
  summarise_sight_cohort(df_prepared %>% filter(scd_status == "SCD"), "All SCD Cases"),
  summarise_sight_cohort(df_prepared %>% filter(scd_status == "No_SCD"), "All No SCD Controls"),
  map_dfr(sort(unique(as.character(df_prepared[[POPULATION_COLUMN_NAME]]))), function(src) {
    summarise_sight_cohort(df_prepared %>% filter(!!sym(POPULATION_COLUMN_NAME) == src), paste0("Cohort: ", src))
  }),
  missing_audit
)

write_csv(table1_sight_rows, file.path(sight_table_dir, "Table1_PartA_SIGHT_Baseline_Characteristics.csv"))

# ------------------------------------------------------------------------------
# 2. Table 2 (Part A): SIGHT Diagnostic Performance with 95% CIs (TRIPOD+AI)
# ------------------------------------------------------------------------------
eval_sight_stratum <- function(df_sub, stratum_name) {
  tp <- sum(df_sub$scd_status == "SCD" & df_sub$.final_pred_class == "SCD", na.rm = TRUE)
  tn <- sum(df_sub$scd_status == "No_SCD" & df_sub$.final_pred_class == "No_SCD", na.rm = TRUE)
  fp <- sum(df_sub$scd_status == "No_SCD" & df_sub$.final_pred_class == "SCD", na.rm = TRUE)
  fn <- sum(df_sub$scd_status == "SCD" & df_sub$.final_pred_class == "No_SCD", na.rm = TRUE)
  n_tot <- tp + tn + fp + fn
  
  f1_val <- if ((tp + fp + fn) > 0 && tp > 0) (2 * tp) / (2 * tp + fp + fn) else NA_real_
  f2_val <- if ((tp + fp + fn) > 0 && tp > 0) (5 * tp) / (5 * tp + 4 * fn + fp) else NA_real_
  
  m_f1 <- tryCatch(f_meas(df_sub, truth = scd_status, estimate = .final_pred_class, estimator = "macro")$.estimate, error = function(e) NA_real_)
  ll_val <- tryCatch(mn_log_loss(df_sub, truth = scd_status, .pred_SCD)$.estimate, error = function(e) NA_real_)
  
  tibble(
    `Evaluation Stratum` = stratum_name,
    `N (Cases / Controls)` = sprintf("%d (%d / %d)", n_tot, tp + fn, tn + fp),
    `TP / FP / FN / TN` = sprintf("%d / %d / %d / %d", tp, fp, fn, tn),
    `Accuracy (95% CI)` = calc_binom_ci(tp + tn, n_tot),
    `Sensitivity (95% CI)` = calc_binom_ci(tp, tp + fn),
    `Specificity (95% CI)` = calc_binom_ci(tn, tn + fp),
    `PPV (95% CI)` = calc_binom_ci(tp, tp + fp),
    `NPV (95% CI)` = calc_binom_ci(tn, tn + fn),
    `F1-Score` = ifelse(is.na(f1_val), "NA", sprintf("%.4f", f1_val)),
    `F2-Score` = ifelse(is.na(f2_val), "NA", sprintf("%.4f", f2_val)),
    `Macro F1` = ifelse(is.na(m_f1), "NA", sprintf("%.4f", m_f1)),
    `AUROC (95% CI)` = calc_auroc_ci(df_sub$scd_status, df_sub$.pred_SCD),
    `Log Loss` = ifelse(is.na(ll_val), "NA", sprintf("%.4f", ll_val))
  )
}

cv_formatted_for_table <- cv_preds %>%
  rename(.final_pred_class = pred_class)

table2_sight_rows <- bind_rows(
  eval_sight_stratum(cv_formatted_for_table, "10-Fold Cross-Validation (Training Set)"),
  eval_sight_stratum(final_results, "Validation Test Set (Overall)"),
  map_dfr(sort(unique(as.character(final_results[[POPULATION_COLUMN_NAME]]))), function(src) {
    eval_sight_stratum(
      final_results %>% filter(!!sym(POPULATION_COLUMN_NAME) == src),
      paste0("Test Subgroup: ", src)
    )
  })
)

write_csv(table2_sight_rows, file.path(sight_table_dir, "Table2_PartA_SIGHT_Diagnostic_Performance.csv"))

# ------------------------------------------------------------------------------
# 3. Supplementary Table S2 (Part A): SIGHT Algorithm Benchmarking & Hyperparameters
# ------------------------------------------------------------------------------
benchmark_sight_table <- all_threshold_perf %>%
  group_by(wflow_id) %>%
  arrange(desc(f2_score), desc(sensitivity)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Algorithm = recode(wflow_id,
                       "recipe_xgboost"  = "XGBoost",
                       "recipe_lightgbm" = "LightGBM",
                       "recipe_rf"       = "Random Forest",
                       "recipe_logreg"   = "Penalized Logistic Regression"
    ),
    `Optimal Threshold` = sprintf("%.2f", threshold),
    `Screening F2-Score` = sprintf("%.4f", f2_score),
    `Screening Sensitivity` = sprintf("%.1f%%", sensitivity * 100),
    `Screening Specificity` = sprintf("%.1f%%", specificity * 100),
    `Screening Accuracy` = sprintf("%.1f%%", accuracy * 100)
  ) %>%
  select(Algorithm, `Optimal Threshold`, `Screening F2-Score`, `Screening Sensitivity`, `Screening Specificity`, `Screening Accuracy`, `.config`)

write_csv(benchmark_sight_table, file.path(sight_table_dir, "Supp_TableS2_PartA_SIGHT_Algorithm_Benchmarking.csv"))
write_csv(best_params, file.path(sight_table_dir, "Supp_TableS2_PartB_SIGHT_Winning_Hyperparameters.csv"))

cat(sprintf("\n✅ All SIGHT Tables exported to '%s/'\n", sight_table_dir))
cat("\n✅ Two-Stage Pipeline Complete.\n")