# ==============================================================================
# SIGHT: Sickle-cell Identification using General Hematology Testing (v2.0)
# ==============================================================================
# Interactive web application for frontline SCD screening from a single Complete
# Blood Count (CBC) using the latest multi-center LightGBM model with exact
# TreeSHAP explanations for both single patients and uploaded batch spreadsheets.
# Visual design matches SHARP while preserving SIGHT's purple/indigo palette.
# ==============================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(tidyverse)
library(tidymodels)
library(bonsai)
library(lightgbm)
library(xgboost)
library(shapviz)
library(vip)
library(bslib)
library(purrr)

# Resolve namespace conflicts
library(conflicted)
conflicts_prefer(shiny::observe)
conflicts_prefer(shiny::dataTableOutput)
conflicts_prefer(shiny::renderDataTable)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
# 1. LOAD LATEST TRAINED SIGHT MODEL
# ==============================================================================
load_sight_model <- function() {
  if (file.exists("scd_final_lightgbm_model.rds")) {
    readRDS("scd_final_lightgbm_model.rds")
  } else if (file.exists("scd_final_model.rds")) {
    readRDS("scd_final_model.rds")
  } else {
    stop("Required SIGHT model file ('scd_final_lightgbm_model.rds' or 'scd_final_model.rds') not found.")
  }
}

final_model <- load_sight_model()

# Clean CBC feature names for clinical display
clean_sight_feature_names <- function(x) {
  dplyr::recode(
    as.character(x),
    "total_hgb"  = "Total Hemoglobin (g/dL)",
    "mcv"        = "MCV (fL)",
    "mchc"       = "MCHC (g/dL)",
    "rdw"        = "RDW (%)",
    "hematocrit" = "Hematocrit (%)",
    "rbc_count"  = "RBC Count (M/μL)",
    .default     = as.character(x)
  )
}

# ==============================================================================
# 2. VALIDATION, SINGLE PREDICTION, SHAP, AND BATCH ANALYSIS FUNCTIONS
# ==============================================================================

#' Validate SIGHT single-patient CBC inputs
validate_sight_inputs <- function(input) {
  errors <- c()
  required_fields <- list(
    "total_hgb"  = "Total Hemoglobin (g/dL)",
    "mcv"        = "MCV (fL)",
    "mchc"       = "MCHC (g/dL)",
    "rdw"        = "RDW (%)",
    "hematocrit" = "Hematocrit (%)",
    "rbc_count"  = "RBC Count (M/μL)"
  )

  for (field in names(required_fields)) {
    val <- input[[field]]
    if (is.null(val) || is.na(val) || val == "") {
      errors <- c(errors, paste("Missing required field:", required_fields[[field]]))
    }
  }

  if (!is.null(input$total_hgb) && !is.na(input$total_hgb)) {
    if (input$total_hgb < 2 || input$total_hgb > 25) {
      errors <- c(errors, "Total Hemoglobin must be between 2 and 25 g/dL")
    }
  }
  if (!is.null(input$mcv) && !is.na(input$mcv)) {
    if (input$mcv < 40 || input$mcv > 140) {
      errors <- c(errors, "MCV must be between 40 and 140 fL")
    }
  }
  if (!is.null(input$mchc) && !is.na(input$mchc)) {
    if (input$mchc < 20 || input$mchc > 45) {
      errors <- c(errors, "MCHC must be between 20 and 45 g/dL")
    }
  }
  if (!is.null(input$rdw) && !is.na(input$rdw)) {
    if (input$rdw < 8 || input$rdw > 45) {
      errors <- c(errors, "RDW must be between 8% and 45%")
    }
  }
  if (!is.null(input$hematocrit) && !is.na(input$hematocrit)) {
    if (input$hematocrit < 6 || input$hematocrit > 75) {
      errors <- c(errors, "Hematocrit must be between 6% and 75%")
    }
  }
  if (!is.null(input$rbc_count) && !is.na(input$rbc_count)) {
    if (input$rbc_count < 0.8 || input$rbc_count > 10) {
      errors <- c(errors, "RBC Count must be between 0.8 and 10 M/μL")
    }
  }

  list(valid = length(errors) == 0, errors = errors)
}

#' Check clinical hematology warnings
check_sight_clinical_warnings <- function(input_data) {
  warnings <- list()
  if (!is.na(input_data$total_hgb) && input_data$total_hgb < 7.0) {
    warnings <- append(warnings, "Severe anemia detected (Total Hgb < 7.0 g/dL) — urgent clinical evaluation advised.")
  }
  if (!is.na(input_data$rdw) && input_data$rdw > 18.0) {
    warnings <- append(warnings, "Markedly elevated RDW (>18%) indicates pronounced anisocytosis, common in sickle cell syndromes and hemolytic anemias.")
  }
  if (!is.na(input_data$mcv) && input_data$mcv < 72.0) {
    warnings <- append(warnings, "Microcytosis (MCV < 72 fL) — consider co-existing S-β thalassemia, α-thalassemia trait, or iron deficiency anemia.")
  }
  if (!is.na(input_data$mcv) && input_data$mcv > 102.0) {
    warnings <- append(warnings, "Macrocytosis (MCV > 102 fL) — may reflect reticulocytosis, hydroxyurea therapy, or folate/B12 deficiency.")
  }
  # Check physiological Hgb * 3 ~ Hematocrit relationship
  if (!is.na(input_data$total_hgb) && !is.na(input_data$hematocrit)) {
    ratio <- input_data$hematocrit / input_data$total_hgb
    if (ratio < 2.3 || ratio > 3.7) {
      warnings <- append(warnings, sprintf("Unusual Hematocrit-to-Hemoglobin ratio (%.1f) — verify CBC analyzer calibration or MCHC abnormality.", ratio))
    }
  }
  list(has_warnings = length(warnings) > 0, warnings = warnings)
}

#' Build a dataframe matching the SIGHT recipe schema
build_sight_recipe_df <- function(df) {
  n <- nrow(df)
  tibble::tibble(
    patient_mrn = if ("patient_id" %in% names(df)) as.character(df$patient_id) else paste0("PT_", seq_len(n)),
    genotype    = "Unknown",
    source      = "Clinical",
    scd_status  = factor(rep("No_SCD", n), levels = c("No_SCD", "SCD")),
    hematocrit  = as.numeric(df$hematocrit),
    mchc        = as.numeric(df$mchc),
    mcv         = as.numeric(df$mcv),
    rbc_count   = as.numeric(df$rbc_count),
    total_hgb   = as.numeric(df$total_hgb),
    rdw         = as.numeric(df$rdw)
  )
}

#' Compute SIGHT prediction and exact TreeSHAP for a single patient
run_sight_single_prediction <- function(input_data, model, threshold = 0.35) {
  single_df <- tibble::tibble(
    patient_id = input_data$patient_id %||% "Single_Patient",
    total_hgb  = as.numeric(input_data$total_hgb),
    mcv        = as.numeric(input_data$mcv),
    mchc       = as.numeric(input_data$mchc),
    rdw        = as.numeric(input_data$rdw),
    hematocrit = as.numeric(input_data$hematocrit),
    rbc_count  = as.numeric(input_data$rbc_count)
  )

  recipe_df <- build_sight_recipe_df(single_df)
  pred_prob <- predict(model, new_data = recipe_df, type = "prob")
  scd_prob  <- as.numeric(pred_prob$.pred_SCD[1])
  no_scd_prob <- as.numeric(pred_prob$.pred_No_SCD[1])

  is_positive <- scd_prob >= threshold
  pred_status <- ifelse(is_positive, "SCD", "No_SCD")
  screening_call <- ifelse(is_positive, "SCD POSITIVE", "SCD NEGATIVE")

  # Extract fitted recipe and engine for exact TreeSHAP
  fitted_rec <- workflows::extract_recipe(model, estimated = TRUE)
  baked_x <- recipes::bake(fitted_rec, new_data = recipe_df, recipes::all_predictors())
  mat_x <- as.matrix(baked_x)
  engine_fit <- workflows::extract_fit_engine(model)

  sv <- shapviz::shapviz(engine_fit, X_pred = mat_x, X = as.data.frame(baked_x))
  if (inherits(sv, "mshapviz") || (is.list(sv) && !inherits(sv, "shapviz"))) {
    sv <- sv[[1]]
  }

  shap_vec <- as.numeric(sv$S[1, ])
  names(shap_vec) <- colnames(sv$S)
  baseline_val <- as.numeric(sv$baseline)

  raw_vals <- c(
    "hematocrit" = single_df$hematocrit,
    "mchc"       = single_df$mchc,
    "mcv"        = single_df$mcv,
    "rbc_count"  = single_df$rbc_count,
    "total_hgb"  = single_df$total_hgb,
    "rdw"        = single_df$rdw
  )

  contrib_df <- tibble::tibble(
    feature = colnames(sv$S),
    feature_label_clean = clean_sight_feature_names(colnames(sv$S)),
    feature_value = as.numeric(raw_vals[colnames(sv$S)]),
    shap_value = shap_vec,
    abs_shap = abs(shap_vec)
  ) %>%
    dplyr::arrange(dplyr::desc(abs_shap)) %>%
    dplyr::mutate(
      rank = dplyr::row_number(),
      direction = ifelse(shap_value >= 0, "Positive", "Negative"),
      clinical_interpretation = purrr::pmap_chr(
        list(feature_label_clean, feature_value, shap_value),
        function(lbl, val, s_val) {
          dir_str <- ifelse(s_val >= 0, "increases", "decreases")
          mag_str <- ifelse(abs(s_val) > 1.0, "strongly", ifelse(abs(s_val) > 0.3, "moderately", "mildly"))
          sprintf("%s = %.1f %s %s predicted SCD risk (SHAP = %+.3f)", lbl, val, mag_str, dir_str, s_val)
        }
      )
    )

  waterfall_plot <- create_sight_waterfall_plot(contrib_df, baseline_val, scd_prob)
  importance_plot <- create_sight_importance_plot(contrib_df)
  warnings <- check_sight_clinical_warnings(single_df)

  list(
    prediction_result = list(
      predicted_status = pred_status,
      screening_call   = screening_call,
      scd_prob         = scd_prob,
      no_scd_prob      = no_scd_prob,
      threshold        = threshold,
      is_positive      = is_positive,
      clinical_warnings = warnings
    ),
    shap_data = list(
      waterfall_plot        = waterfall_plot,
      importance_plot       = importance_plot,
      feature_contributions = contrib_df,
      baseline              = baseline_val,
      shap_values           = shap_vec
    ),
    input_data = single_df
  )
}

#' Create horizontal SHAP waterfall plot for SIGHT (matching SHARP style)
create_sight_waterfall_plot <- function(feature_contributions, baseline_prob = 0.0, scd_prob = NA) {
  plot_data <- feature_contributions %>%
    dplyr::arrange(abs_shap)

  cumulative_values <- cumsum(c(baseline_prob, plot_data$shap_value))
  plot_data$start_pos <- cumulative_values[-length(cumulative_values)]
  plot_data$end_pos <- cumulative_values[-1]

  plot_data <- plot_data %>%
    dplyr::mutate(
      feature_label = paste0(feature_label_clean, " = ", round(feature_value, 1)),
      y_idx = dplyr::row_number()
    )

  final_log_odds <- baseline_prob + sum(plot_data$shap_value)

  ggplot(plot_data) +
    geom_rect(aes(
      xmin = pmin(start_pos, end_pos),
      xmax = pmax(start_pos, end_pos),
      ymin = y_idx - 0.35,
      ymax = y_idx + 0.35,
      fill = direction
    ), color = "white", linewidth = 0.5) +
    geom_text(aes(
      x = (start_pos + end_pos) / 2,
      y = y_idx,
      label = sprintf("%+.3f", shap_value)
    ), color = "#1B1919", fontface = "bold", size = 3.3) +
    geom_vline(xintercept = baseline_prob, linetype = "dashed", color = "gray40", alpha = 0.8) +
    geom_vline(xintercept = final_log_odds, linetype = "solid", color = "#593196", alpha = 0.9) +
    scale_fill_manual(
      values = c("Positive" = "#FFA500", "Negative" = "#593196"),
      labels = c("Positive" = "Increases SCD Risk (+)", "Negative" = "Decreases SCD Risk (-)")
    ) +
    scale_y_continuous(
      breaks = plot_data$y_idx,
      labels = plot_data$feature_label,
      expand = expansion(mult = c(0.10, 0.14))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
    labs(
      title = "SIGHT SHAP Waterfall Plot (CBC Feature Contributions)",
      subtitle = sprintf(
        "Base log-odds E[f(X)] = %.2f  ->  Final log-odds f(x) = %.2f  (Predicted SCD Risk = %.1f%%)",
        baseline_prob, final_log_odds, scd_prob * 100
      ),
      x = "Model Output Contribution (log-odds of SCD)",
      y = NULL,
      fill = "Contribution"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#' Create patient-specific SHAP feature importance plot for SIGHT
create_sight_importance_plot <- function(feature_contributions) {
  ggplot(feature_contributions, aes(x = reorder(feature_label_clean, abs_shap), y = abs_shap, fill = direction)) +
    geom_col(width = 0.68, color = "white", linewidth = 0.4) +
    coord_flip() +
    scale_fill_manual(
      values = c("Positive" = "#d9534f", "Negative" = "#593196"),
      labels = c("Positive" = "Increases SCD Risk", "Negative" = "Decreases SCD Risk")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(
      title = "CBC Feature Importance (|SHAP Value|)",
      subtitle = "Magnitude of each CBC parameter's impact on this patient's SCD screening score",
      x = NULL,
      y = "|SHAP Value| (log-odds impact)",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    )
}

#' Summarize SIGHT SHAP insights for the Explanation tab
summarize_sight_shap_insights <- function(contrib_df, pred_res) {
  insights <- c()
  top_f <- contrib_df[1, ]
  insights <- c(insights, paste0("🔍 Primary hematological driver: ", top_f$clinical_interpretation))

  pos_df <- contrib_df %>% dplyr::filter(shap_value > 0)
  neg_df <- contrib_df %>% dplyr::filter(shap_value < 0)

  if (nrow(pos_df) > 0) {
    insights <- c(insights, paste0("⬆️ Risk-increasing CBC features: ", paste(pos_df$feature_label_clean, collapse = ", ")))
  }
  if (nrow(neg_df) > 0) {
    insights <- c(insights, paste0("⬇️ Risk-decreasing CBC features: ", paste(neg_df$feature_label_clean, collapse = ", ")))
  }

  if (pred_res$is_positive) {
    insights <- c(insights, sprintf(
      "⚠️ Predicted SCD probability (%.1f%%) exceeds screening threshold (%.0f%%) — confirmatory hemoglobin electrophoresis/HPLC (SHARP) is recommended.",
      pred_res$scd_prob * 100, pred_res$threshold * 100
    ))
  } else {
    insights <- c(insights, sprintf(
      "✅ Predicted SCD probability (%.1f%%) is below screening threshold (%.0f%%) — hematological profile is unlikely to represent Sickle Cell Disease.",
      pred_res$scd_prob * 100, pred_res$threshold * 100
    ))
  }
  insights
}

# ==============================================================================
# 3. BATCH SPREADSHEET ANALYSIS MODULE FOR SIGHT (WITH FULL SHAP COLUMNS)
# ==============================================================================

SIGHT_BATCH_COLUMN_MAP <- list(
  patient_id = c("patient_id", "patient_mrn", "mrn", "id", "sample_id", "subject_id"),
  total_hgb  = c("total_hgb", "hgb", "hemoglobin", "hemoglobin-blood", "hb", "total hemoglobin", "hgb (g/dl)"),
  mcv        = c("mcv", "mean cell volume", "mean corpuscular volume", "mcv-mean cell volume fl", "mcv (fl)"),
  mchc       = c("mchc", "mean cell hemoglobin concentration", "mean corpuscular hemoglobin concentration", "mchc (g/dl)"),
  rdw        = c("rdw", "red cell distribution width", "rdw (%)", "rdw-cv"),
  hematocrit = c("hematocrit", "hct", "hematocrit-blood", "hct (%)"),
  rbc_count  = c("rbc_count", "rbc", "red blood cell count", "rbc count-blood", "rbc (m/ul)", "rbc (x10^12/l)")
)

#' Generate template dataframe for SIGHT batch upload
get_sight_batch_template <- function() {
  tibble::tibble(
    patient_id = c("PT_001_Classic_SS", "PT_002_Severe_SCD", "PT_003_SC_Disease", "PT_004_Sbeta0", "PT_005_Healthy_Control", "PT_006_Normal_Adult"),
    total_hgb  = c(7.6,  6.4,  10.8, 8.2,  14.2, 13.5),
    mcv        = c(91.0, 88.5, 78.0, 67.5, 89.0, 87.5),
    mchc       = c(34.5, 33.8, 35.2, 32.0, 33.8, 34.0),
    rdw        = c(22.4, 24.1, 16.8, 20.5, 12.6, 13.1),
    hematocrit = c(22.5, 19.0, 31.5, 25.0, 42.5, 40.2),
    rbc_count  = c(2.48, 2.15, 4.05, 3.70, 4.78, 4.60)
  )
}

#' Standardize SIGHT batch columns
standardize_sight_batch_columns <- function(df) {
  current_names <- names(df)
  names_clean <- trimws(tolower(gsub("[_\\.]+", " ", current_names)))
  names_exact <- trimws(tolower(current_names))

  std_df <- df
  mapped_info <- c()

  for (std_name in names(SIGHT_BATCH_COLUMN_MAP)) {
    aliases <- tolower(SIGHT_BATCH_COLUMN_MAP[[std_name]])
    aliases_clean <- gsub("[_\\.]+", " ", aliases)
    idx <- which(names_exact %in% aliases | names_clean %in% aliases_clean)
    if (length(idx) > 0) {
      orig_col <- current_names[idx[1]]
      std_df[[std_name]] <- df[[orig_col]]
      mapped_info <- c(mapped_info, sprintf("%s -> %s", orig_col, std_name))
    }
  }

  required_cols <- c("total_hgb", "mcv", "mchc", "rdw", "hematocrit", "rbc_count")
  missing_req <- setdiff(required_cols, names(std_df))

  list(
    data = std_df,
    mapped = mapped_info,
    missing_required = missing_req
  )
}

#' Run SIGHT batch prediction + exact LightGBM TreeSHAP on an uploaded dataframe
run_sight_batch_analysis <- function(raw_df, model, threshold = 0.35, progress_callback = NULL) {
  if (is.null(raw_df) || nrow(raw_df) == 0) {
    return(list(results_df = NULL, summary_stats = NULL, error = "Uploaded spreadsheet is empty."))
  }

  if (!is.null(progress_callback)) progress_callback(0.15, "Validating and standardizing CBC columns...")
  std_res <- standardize_sight_batch_columns(raw_df)
  if (length(std_res$missing_required) > 0) {
    return(list(
      results_df = NULL,
      summary_stats = NULL,
      error = paste0(
        "Missing required CBC columns: ",
        paste(std_res$missing_required, collapse = ", "),
        ". Please download the Template CSV to see the exact required column format."
      )
    ))
  }

  df <- std_res$data
  n_rows <- nrow(df)

  cbc_cols <- c("total_hgb", "mcv", "mchc", "rdw", "hematocrit", "rbc_count")
  for (col in cbc_cols) {
    df[[col]] <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", as.character(df[[col]]))))
  }

  if (!("patient_id" %in% names(df))) {
    df$patient_id <- paste0("Row_", seq_len(n_rows))
  }

  if (!is.null(progress_callback)) progress_callback(0.40, "Running LightGBM SCD screening predictions...")
  recipe_df <- build_sight_recipe_df(df)
  probs <- predict(model, new_data = recipe_df, type = "prob")

  scd_prob <- round(as.numeric(probs$.pred_SCD), 4)
  no_scd_prob <- round(as.numeric(probs$.pred_No_SCD), 4)
  is_pos <- !is.na(scd_prob) & (scd_prob >= threshold)
  pred_status <- ifelse(is.na(scd_prob), "Incomplete_CBC", ifelse(is_pos, "SCD", "No_SCD"))
  screening_call <- ifelse(is.na(scd_prob), "INCOMPLETE DATA", ifelse(is_pos, "SCD POSITIVE", "SCD NEGATIVE"))
  rec_action <- ifelse(
    is.na(scd_prob),
    "Verify missing CBC parameters",
    ifelse(is_pos, "Confirmatory Electrophoresis / HPLC (Run SHARP)", "Low Risk — Routine Clinical Follow-Up")
  )

  if (!is.null(progress_callback)) progress_callback(0.70, "Computing exact LightGBM TreeSHAP values for all rows...")
  fitted_rec <- workflows::extract_recipe(model, estimated = TRUE)
  baked_x <- recipes::bake(fitted_rec, new_data = recipe_df, recipes::all_predictors())
  mat_x <- as.matrix(baked_x)
  engine_fit <- workflows::extract_fit_engine(model)

  sv <- shapviz::shapviz(engine_fit, X_pred = mat_x, X = as.data.frame(baked_x))
  if (inherits(sv, "mshapviz") || (is.list(sv) && !inherits(sv, "shapviz"))) {
    sv <- sv[[1]]
  }

  shap_mat <- as.matrix(sv$S)
  baseline_val <- as.numeric(sv$baseline)
  clean_col_labels <- clean_sight_feature_names(colnames(shap_mat))

  top_pos_vec <- character(n_rows)
  top_neg_vec <- character(n_rows)
  summary_vec <- character(n_rows)

  for (i in seq_len(n_rows)) {
    row_s <- shap_mat[i, ]
    if (all(is.na(row_s))) {
      top_pos_vec[i] <- "NA"
      top_neg_vec[i] <- "NA"
      summary_vec[i] <- "Missing CBC values"
      next
    }
    pos_idx <- which(row_s > 0)
    if (length(pos_idx) > 0) {
      best_p <- pos_idx[which.max(row_s[pos_idx])]
      top_pos_vec[i] <- sprintf("%s (%+.3f)", clean_col_labels[best_p], row_s[best_p])
    } else {
      top_pos_vec[i] <- "None"
    }

    neg_idx <- which(row_s < 0)
    if (length(neg_idx) > 0) {
      best_n <- neg_idx[which.min(row_s[neg_idx])]
      top_neg_vec[i] <- sprintf("%s (%+.3f)", clean_col_labels[best_n], row_s[best_n])
    } else {
      top_neg_vec[i] <- "None"
    }

    top1 <- which.max(abs(row_s))
    dir_txt <- ifelse(row_s[top1] >= 0, "increases SCD risk", "decreases SCD risk")
    summary_vec[i] <- sprintf("Primary driver: %s (%+.3f, %s)", clean_col_labels[top1], row_s[top1], dir_txt)
  }

  shap_df <- as.data.frame(round(shap_mat, 4))
  get_shap_col <- function(nm) if (nm %in% names(shap_df)) shap_df[[nm]] else rep(NA_real_, n_rows)

  out_cols <- tibble::tibble(
    predicted_status          = pred_status,
    screening_call            = screening_call,
    scd_probability           = scd_prob,
    likelihood_of_SCD_percent = round(scd_prob * 100, 2),
    no_scd_probability        = no_scd_prob,
    threshold_used            = threshold,
    recommended_action        = rec_action,
    shap_baseline_log_odds    = round(baseline_val, 4),
    shap_prediction_log_odds  = round(baseline_val + rowSums(shap_mat, na.rm = TRUE), 4),
    shap_total_hgb            = get_shap_col("total_hgb"),
    shap_mcv                  = get_shap_col("mcv"),
    shap_mchc                 = get_shap_col("mchc"),
    shap_rdw                  = get_shap_col("rdw"),
    shap_hematocrit           = get_shap_col("hematocrit"),
    shap_rbc_count            = get_shap_col("rbc_count"),
    top_positive_shap_feature = top_pos_vec,
    top_negative_shap_feature = top_neg_vec,
    shap_clinical_summary     = summary_vec
  )

  results_df <- dplyr::bind_cols(raw_df, out_cols)

  summary_stats <- list(
    total_rows     = n_rows,
    scd_positive   = sum(pred_status == "SCD", na.rm = TRUE),
    scd_negative   = sum(pred_status == "No_SCD", na.rm = TRUE),
    positive_rate  = round(mean(pred_status == "SCD", na.rm = TRUE) * 100, 1),
    mean_scd_prob  = round(mean(scd_prob, na.rm = TRUE) * 100, 1),
    threshold_used = threshold
  )

  if (!is.null(progress_callback)) progress_callback(1.0, "Batch analysis complete!")

  list(
    results_df = results_df,
    summary_stats = summary_stats,
    error = NULL
  )
}

# ==============================================================================
# 4. UI RENDERING HELPERS FOR SIGHT (MATCHING SHARP STRUCTURE)
# ==============================================================================

render_sight_prediction_results <- function(pred_res, input_data) {
  if (is.null(pred_res)) {
    return(div(
      class = "alert alert-info",
      icon("circle-info"),
      h4("Ready for SIGHT Screening"),
      p("Enter a patient's 6 routine Complete Blood Count (CBC) parameters in the sidebar and click 'Predict SCD Status' to begin analysis.")
    ))
  }

  card_border <- if (pred_res$is_positive) "border-danger" else "border-success"
  badge_class <- if (pred_res$is_positive) "bg-danger" else "bg-success"
  badge_text  <- if (pred_res$is_positive) "Confirmatory Genotyping Recommended (SHARP)" else "Low Risk — Automated Screen Negative"
  heading_col <- if (pred_res$is_positive) "text-danger" else "text-success"

  prob_tbl <- tibble::tibble(
    `Screening Class` = c("SCD (Sickle Cell Disease)", "No_SCD (Unlikely SCD / Control)"),
    `Predicted Probability` = c(
      paste0(round(pred_res$scd_prob * 100, 2), "%"),
      paste0(round(pred_res$no_scd_prob * 100, 2), "%")
    )
  )

  ref_tbl <- tibble::tibble(
    `CBC Parameter` = c("Total Hemoglobin", "Hematocrit", "MCV", "MCHC", "RDW", "RBC Count"),
    `Patient Value` = c(
      paste0(input_data$total_hgb, " g/dL"),
      paste0(input_data$hematocrit, " %"),
      paste0(input_data$mcv, " fL"),
      paste0(input_data$mchc, " g/dL"),
      paste0(input_data$rdw, " %"),
      paste0(input_data$rbc_count, " M/μL")
    ),
    `Typical Normal Range` = c("12.0–16.0 g/dL", "36.0–48.0 %", "80.0–100.0 fL", "32.0–36.0 g/dL", "11.5–14.5 %", "4.0–5.5 M/μL"),
    `Typical SCD Pattern`  = c("6.0–10.5 g/dL", "18.0–30.0 %", "70.0–98.0 fL", "32.0–37.0 g/dL", "16.0–26.0 %", "2.0–4.0 M/μL")
  )

  div(
    class = "row",
    div(
      class = "col-12",
      div(
        class = paste("card", card_border),
        div(
          class = "card-header d-flex justify-content-between align-items-center",
          h4("Screening Result", class = "mb-0"),
          span(class = paste("badge", badge_class), badge_text)
        ),
        div(
          class = "card-body",
          div(
            class = "row mb-3",
            div(
              class = "col-md-6",
              h2(
                if (pred_res$is_positive) "SCD POSITIVE (High SCD Likelihood)" else "SCD NEGATIVE (Low SCD Likelihood)",
                class = paste(heading_col, "mb-2 fw-bold")
              ),
              p(class = "lead mb-1", paste0("Predicted SCD Probability: ", round(pred_res$scd_prob * 100, 1), "%")),
              p(class = "text-muted mb-0", paste0("Decision Threshold: ", round(pred_res$threshold * 100, 1), "% (LightGBM Multi-Center Classifier)"))
            ),
            div(
              class = "col-md-6",
              h5("Clinical Interpretation & Next Steps:"),
              p(
                if (pred_res$is_positive) {
                  "This patient's CBC pattern strongly resembles Sickle Cell Disease (e.g., hemolytic anemia with elevated RDW and reduced RBC/Hemoglobin). Proceed to Stage 2 confirmatory testing via Hemoglobin Electrophoresis / HPLC and run SHARP for definitive 5-class genotyping."
                } else {
                  "This patient's CBC profile does not exhibit the characteristic hematological signature of Sickle Cell Disease. Note: Silent carrier states (such as Sickle Cell Trait, HbAS) typically present with normal CBC indices and require electrophoresis if carrier screening is clinically indicated."
                }
              )
            )
          ),

          if (pred_res$is_positive) {
            div(
              class = "alert alert-danger mb-3",
              icon("triangle-exclamation"),
              strong(" Action Recommended — Stage 2 Confirmatory Testing (SHARP):"),
              br(),
              "Order quantitative Hemoglobin Electrophoresis or HPLC (HbS%, HbA%, HbC%, HbF%, HbA2%) and evaluate with SHARP to distinguish HbSS, HbSC, HbSβ⁰, and HbSβ⁺."
            )
          },

          if (!is.null(pred_res$clinical_warnings) && pred_res$clinical_warnings$has_warnings) {
            div(
              class = "alert alert-warning mb-3",
              icon("exclamation-triangle"),
              strong(" Hematological Considerations:"),
              tags$ul(
                class = "mb-0 mt-1",
                lapply(pred_res$clinical_warnings$warnings, function(w) tags$li(w))
              )
            )
          }
        )
      )
    ),

    div(
      class = "col-12 mt-3",
      div(
        class = "card",
        div(class = "card-header", h5("Predicted Class Probabilities", class = "mb-0")),
        div(
          class = "card-body p-2",
          DT::datatable(
            prob_tbl,
            options = list(dom = "t", paging = FALSE, searching = FALSE, ordering = FALSE),
            rownames = FALSE
          )
        )
      )
    ),

    div(
      class = "col-12 mt-3",
      div(
        class = "card",
        div(class = "card-header", h5("Patient CBC vs. Clinical Reference Ranges", class = "mb-0")),
        div(
          class = "card-body p-2",
          DT::datatable(
            ref_tbl,
            options = list(dom = "t", paging = FALSE, searching = FALSE, ordering = FALSE),
            rownames = FALSE
          )
        )
      )
    )
  )
}

render_sight_shap_explanation <- function(shap_data, pred_res) {
  if (is.null(shap_data)) {
    return(div(
      class = "alert alert-info",
      icon("circle-info"),
      h4("SHAP Explanation"),
      p("Patient-specific TreeSHAP feature contributions and waterfall visualizations will appear here after you click 'Predict SCD Status'.")
    ))
  }

  div(
    div(
      class = "mb-4",
      h4("CBC Feature Contributions (TreeSHAP Waterfall)"),
      plotOutput("shap_waterfall", height = "380px")
    ),
    div(
      class = "mb-4",
      h4("Key Hematological Insights"),
      div(
        class = "alert alert-light",
        tags$ul(
          class = "mb-0",
          lapply(summarize_sight_shap_insights(shap_data$feature_contributions, pred_res), function(ins) {
            tags$li(HTML(ins))
          })
        )
      )
    ),
    div(
      h4("Detailed CBC Feature Contributions"),
      DT::dataTableOutput("shap_contributions_table")
    )
  )
}

render_sight_clinical_guidelines <- function() {
  div(
    h3("SIGHT Clinical Screening Guidelines"),

    div(
      class = "card guideline-card mb-3",
      div(class = "card-header", h5("Characteristic CBC Signatures in Sickle Cell Disease")),
      div(
        class = "card-body",
        tags$ul(
          tags$li(strong("RDW (%): "), "Typically markedly elevated (>16–25%) reflecting chronic hemolysis, reticulocytosis, and poikilocytosis/sickled forms."),
          tags$li(strong("Total Hemoglobin & Hematocrit: "), "Moderately to severely reduced (Hgb 6–10 g/dL, HCT 18–30%) in HbSS and HbSβ⁰; mildly reduced (9–12 g/dL) in HbSC and HbSβ⁺."),
          tags$li(strong("RBC Count: "), "Reduced in proportion to hemolytic anemia (typically 2.0–3.5 M/μL in HbSS), whereas thalassemia trait often shows elevated RBC count with microcytosis."),
          tags$li(strong("MCV (fL): "), "Normocytic (80–95 fL) or mildly macrocytic in classic HbSS; microcytic (<75 fL) in HbSβ⁰, HbSβ⁺, or co-inherited α-thalassemia."),
          tags$li(strong("MCHC (g/dL): "), "Often normal to upper-normal due to cellular dehydration (especially in HbSC disease).")
        )
      )
    ),

    div(
      class = "alert alert-info",
      h5("Important Screening Considerations"),
      tags$ul(
        class = "mb-0",
        tags$li(strong("Sickle Cell Trait (HbAS): "), "Asymptomatic carriers (HbAS) generally have completely normal CBC parameters and are classified as 'No_SCD' by SIGHT. Carrier detection requires hemoglobin electrophoresis/HPLC."),
        tags$li(strong("Threshold Selection: "), "The default F2-optimized threshold (0.35) balances high sensitivity (>94%) and specificity (>97%). For ultra-high-sensitivity rule-out triage, the 0.10 threshold provides >97.8% sensitivity.")
      )
    )
  )
}

render_sight_about_content <- function() {
  div(
    h3("About SIGHT"),

    div(
      class = "card mb-3",
      div(class = "card-header", h5("Overview")),
      div(
        class = "card-body",
        p("SIGHT (Sickle-cell Identification using General Hematology Testing) is a machine learning clinical decision support tool that screens for Sickle Cell Disease (SCD) from a single routine Complete Blood Count (CBC)."),
        p("Designed to expand rapid screening access in both acute emergency care and resource-limited global health settings (LMICs), SIGHT requires no specialized equipment beyond a standard hematology analyzer.")
      )
    ),

    div(
      class = "card mb-3",
      div(class = "card-header", h5("Latest LightGBM Model Information")),
      div(
        class = "card-body",
        tags$ul(
          tags$li(strong("Algorithm: "), "Gradient Boosted Decision Trees (LightGBM via tidymodels & bonsai)"),
          tags$li(strong("Multi-Center Training Cohort: "), "~9,900 cleaned CBC profiles across US (Grady Memorial Hospital & NHANES), Ghanaian, and Senegalese cohorts"),
          tags$li(strong("Predictors (6 CBC Features): "), "Total Hemoglobin, Hematocrit, RBC Count, MCV, MCHC, and RDW"),
          tags$li(strong("Validation Performance: "), "AUROC > 0.991; 96.6% overall accuracy (94.1% sensitivity, 97.2% specificity) at the F2-optimized 0.35 cutoff"),
          tags$li(strong("Explainability: "), "Exact TreeSHAP (shapviz) log-odds feature attributions for both single-patient and batch spreadsheet analyses")
        )
      )
    ),

    div(
      class = "card mb-3",
      div(class = "card-header", h5("Important Disclaimers")),
      div(
        class = "card-body",
        div(
          class = "alert alert-warning mb-0",
          tags$ul(
            class = "mb-0",
            tags$li("This tool is intended for research, educational, and clinical decision support purposes."),
            tags$li("SIGHT is a screening tool and does not replace confirmatory hemoglobin electrophoresis, HPLC, or genetic testing."),
            tags$li("Other severe hemolytic or microcytic anemias can share CBC features with SCD and warrant confirmatory testing.")
          )
        )
      )
    )
  )
}

generate_sight_report <- function(pred_res, shap_data, input_data, file) {
  if (is.null(pred_res)) {
    writeLines("<html><body><h2>No prediction has been run yet.</h2></body></html>", file)
    return()
  }

  contrib_rows <- ""
  if (!is.null(shap_data) && !is.null(shap_data$feature_contributions)) {
    df_c <- shap_data$feature_contributions
    contrib_rows <- paste0(
      sapply(seq_len(nrow(df_c)), function(i) {
        sprintf(
          "<tr><td>%s</td><td>%.2f</td><td>%+.4f</td><td>%s</td><td>%s</td></tr>",
          df_c$feature_label_clean[i], df_c$feature_value[i], df_c$shap_value[i],
          df_c$direction[i], df_c$clinical_interpretation[i]
        )
      }),
      collapse = "\n"
    )
  }

  html_content <- paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'><title>SIGHT Screening Report</title>",
    "<style>body{font-family:Arial,sans-serif;margin:40px;color:#2c3e50;}",
    ".header{background:#593196;color:white;padding:20px;border-radius:8px;}",
    ".result{background:#f8f9fa;padding:20px;border-left:5px solid ", ifelse(pred_res$is_positive, "#d9534f", "#2e7d32"), ";border-radius:5px;margin:15px 0;}",
    "table{border-collapse:collapse;width:100%;margin-top:10px;}",
    "th,td{border:1px solid #ddd;padding:10px;text-align:left;}",
    "th{background-color:#ede7f6;color:#593196;}</style></head><body>",
    "<div class='header'><h1>SIGHT Frontline SCD Screening Report</h1>",
    "<p>Sickle-cell Identification using General Hematology Testing (LightGBM + TreeSHAP)</p></div>",
    "<p><strong>Generated on:</strong> ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</p>",
    "<p><strong>Patient ID:</strong> ", input_data$patient_id, "</p>",
    "<h2>Input CBC Parameters</h2>",
    "<table><tr><th>Parameter</th><th>Value</th></tr>",
    "<tr><td>Total Hemoglobin (g/dL)</td><td>", input_data$total_hgb, "</td></tr>",
    "<tr><td>Hematocrit (%)</td><td>", input_data$hematocrit, "</td></tr>",
    "<tr><td>RBC Count (M/uL)</td><td>", input_data$rbc_count, "</td></tr>",
    "<tr><td>MCV (fL)</td><td>", input_data$mcv, "</td></tr>",
    "<tr><td>MCHC (g/dL)</td><td>", input_data$mchc, "</td></tr>",
    "<tr><td>RDW (%)</td><td>", input_data$rdw, "</td></tr></table>",
    "<h2>Screening Result</h2>",
    "<div class='result'>",
    "<h3>", pred_res$screening_call, " (", pred_res$predicted_status, ")</h3>",
    "<p><strong>Predicted SCD Probability:</strong> ", round(pred_res$scd_prob * 100, 2), "%</p>",
    "<p><strong>Classification Threshold:</strong> ", round(pred_res$threshold * 100, 1), "%</p>",
    "</div>",
    "<h2>TreeSHAP Feature Contributions</h2>",
    "<table><tr><th>CBC Feature</th><th>Patient Value</th><th>SHAP (log-odds)</th><th>Direction</th><th>Clinical Interpretation</th></tr>",
    contrib_rows,
    "</table></body></html>"
  )
  writeLines(html_content, file)
}

# ==============================================================================
# 5. USER INTERFACE (SHARP LAYOUT + SIGHT PURPLE THEME)
# ==============================================================================

ui <- page_sidebar(
  title = "SIGHT 👁️🩸",
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#6f42c1",
    secondary = "#95a5a6"
  ),

  tags$head(
    tags$style(HTML("
      /* SIGHT signature deep purple navbar matching SHARP's header design */
      .navbar {
        background-color: #593196 !important;
      }
      .navbar-brand {
        color: white !important;
        font-size: 2.2rem !important;
        font-weight: 600 !important;
      }
    ")),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),

  sidebar = sidebar(
    width = 400,

    # Mode Selector Card
    card(
      card_header("Analysis Mode", class = "bg-sight-primary"),
      radioButtons(
        "analysis_mode",
        label = NULL,
        choices = c(
          "Single Patient Screening" = "single",
          "Batch Spreadsheet Analysis (CSV / Excel)" = "batch"
        ),
        selected = "single"
      )
    ),

    # ---------------- SINGLE PATIENT CONTROLS ----------------
    conditionalPanel(
      condition = "input.analysis_mode == 'single'",

      card(
        card_header("Patient Information"),
        textInput("patient_id", "Patient ID (optional)",
                  placeholder = "Enter patient identifier")
      ),

      card(
        card_header("CBC Hematology Results", class = "bg-sight-primary"),
        helpText("All 6 routine CBC parameters are required (*)"),

        numericInput("total_hgb", "Total Hemoglobin (g/dL)*",
                     value = 14.0, min = 2, max = 25, step = 0.1),
        numericInput("hematocrit", "Hematocrit (%)*",
                     value = 42.0, min = 6, max = 75, step = 0.1),
        numericInput("rbc_count", "RBC Count (M/μL or ×10¹²/L)*",
                     value = 4.8, min = 0.8, max = 10, step = 0.1),
        numericInput("mcv", "MCV (fL)*",
                     value = 90.0, min = 40, max = 140, step = 0.1),
        numericInput("mchc", "MCHC (g/dL)*",
                     value = 34.0, min = 20, max = 45, step = 0.1),
        numericInput("rdw", "RDW (%)*",
                     value = 13.0, min = 8, max = 45, step = 0.1),

        uiOutput("cbc_validation")
      ),

      card(
        card_header("Screening Threshold Settings"),
        selectInput(
          "threshold_preset",
          "Operating Cutoff",
          choices = c(
            "0.35 — F2-Optimized Default (96.6% Acc, 94.1% Sens)" = "0.35",
            "0.25 — High Sensitivity Balanced (96.0% Acc, 95.4% Sens)" = "0.25",
            "0.10 — Maximum Sensitivity Triage (97.9% Sens)" = "0.10",
            "0.50 — Standard Probability Cutoff" = "0.50"
          ),
          selected = "0.35"
        )
      ),

      div(
        class = "d-grid gap-2",
        actionButton("predict_btn", "Predict SCD Status",
                     class = "btn-primary btn-lg", icon = icon("calculator")),
        actionButton("reset_btn", "Reset Form",
                     class = "btn-secondary"),
        br(),
        downloadButton("download_report", "Download Report",
                       class = "btn-outline-sight")
      )
    ),

    # ---------------- BATCH ANALYSIS CONTROLS ----------------
    conditionalPanel(
      condition = "input.analysis_mode == 'batch'",

      card(
        card_header("1. Format Template", class = "bg-sight-primary"),
        p(class = "small text-muted mb-2",
          "Download the formatted template CSV containing the 6 required CBC headers and 6 sample patient profiles."),
        downloadButton("download_batch_template", "Download Template CSV",
                       class = "btn-outline-sight w-100")
      ),

      card(
        card_header("2. Upload CBC Spreadsheet", class = "bg-sight-primary"),
        fileInput(
          "batch_file",
          "Select CSV or Excel (.xlsx / .xls)",
          accept = c(".csv", ".xlsx", ".xls")
        ),
        selectInput(
          "batch_threshold_preset",
          "Batch Screening Cutoff",
          choices = c(
            "0.35 — F2-Optimized Default" = "0.35",
            "0.25 — High Sensitivity Balanced" = "0.25",
            "0.10 — Maximum Sensitivity Triage" = "0.10",
            "0.50 — Standard Probability Cutoff" = "0.50"
          ),
          selected = "0.35"
        ),
        helpText("Required columns: total_hgb, mcv, mchc, rdw, hematocrit, rbc_count. Optional: patient_id."),
        div(
          class = "d-grid gap-2 mt-2",
          actionButton("run_batch_btn", "Run Batch Analysis + SHAP",
                       class = "btn-primary btn-lg", icon = icon("gears")),
          actionButton("load_demo_batch_btn", "Load Demo Batch (6 Cases)",
                       class = "btn-outline-secondary btn-sm")
        )
      ),

      card(
        card_header("3. Export Results + SHAP"),
        p(class = "small text-muted mb-2",
          "Output spreadsheet includes all original columns, SCD probabilities, screening calls, and exact LightGBM TreeSHAP values for all 6 CBC features."),
        downloadButton("download_batch_results", "Download Results CSV (with SHAP)",
                       class = "btn-success w-100")
      )
    )
  ),

  # Main 2-column content area matching SHARP
  layout_columns(

    # Left Column: Single Prediction Results OR Batch Summary & Format Guide
    card(
      full_screen = TRUE,
      card_header(
        uiOutput("left_card_title_ui"),
        class = "bg-sight-primary"
      ),
      uiOutput("left_panel_content_ui")
    ),

    # Right Column: Explanation, Batch Results & Reference Tabs
    card(
      full_screen = TRUE,
      card_header("Explanation, Batch Results & Clinical Reference", class = "bg-sight-primary"),

      navset_card_tab(
        id = "sight_main_tabs",
        nav_panel(
          "Explanation",
          value = "tab_explanation",
          uiOutput("shap_explanation")
        ),
        nav_panel(
          "Batch Analysis",
          value = "tab_batch",
          uiOutput("batch_results_right_ui")
        ),
        nav_panel(
          "Clinical Guidelines",
          value = "tab_guidelines",
          uiOutput("clinical_guidelines")
        ),
        nav_panel(
          "Feature Importance",
          value = "tab_importance",
          plotlyOutput("feature_importance_plot")
        ),
        nav_panel(
          "About",
          value = "tab_about",
          uiOutput("about_content")
        )
      )
    ),

    col_widths = c(6, 6)
  )
)

# ==============================================================================
# 6. SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {

  values <- reactiveValues(
    prediction_result = NULL,
    shap_data = NULL,
    input_data = NULL,
    validation_status = NULL,
    batch_results_df = NULL,
    batch_summary = NULL,
    batch_error = NULL
  )

  # Automatically switch tabs when toggling analysis_mode
  observeEvent(input$analysis_mode, {
    if (identical(input$analysis_mode, "batch")) {
      nav_select("sight_main_tabs", selected = "tab_batch", session = session)
    } else if (identical(input$analysis_mode, "single")) {
      nav_select("sight_main_tabs", selected = "tab_explanation", session = session)
    }
  }, ignoreInit = TRUE)

  output$left_card_title_ui <- renderUI({
    if (identical(input$analysis_mode, "batch")) {
      span(icon("table"), " Batch Screening Summary & Spreadsheet Specification")
    } else {
      span(icon("stethoscope"), " Prediction Results")
    }
  })

  # Real-time CBC input validation
  observe({
    values$validation_status <- validate_sight_inputs(input)
  })

  output$cbc_validation <- renderUI({
    if (!is.null(values$validation_status) && !values$validation_status$valid) {
      div(
        class = "alert alert-warning mt-2",
        HTML(paste("⚠️", paste(values$validation_status$errors, collapse = "<br>")))
      )
    }
  })

  # Single-patient prediction
  observeEvent(input$predict_btn, {
    validation <- validate_sight_inputs(input)
    if (!validation$valid) {
      showNotification(
        paste("Please fix the following errors:", paste(validation$errors, collapse = ", ")),
        type = "error",
        duration = 5
      )
      return()
    }

    withProgress(message = "Calculating SCD risk & TreeSHAP...", value = 0, {
      incProgress(0.2, detail = "Preparing CBC input data")
      input_list <- list(
        patient_id = if (!is.null(input$patient_id) && nzchar(input$patient_id)) input$patient_id else "Single_Patient",
        total_hgb  = input$total_hgb,
        mcv        = input$mcv,
        mchc       = input$mchc,
        rdw        = input$rdw,
        hematocrit = input$hematocrit,
        rbc_count  = input$rbc_count
      )
      th_val <- as.numeric(input$threshold_preset %||% "0.35")

      incProgress(0.6, detail = "Running LightGBM model and TreeSHAP")
      res <- run_sight_single_prediction(input_list, final_model, threshold = th_val)

      values$prediction_result <- res$prediction_result
      values$shap_data         <- res$shap_data
      values$input_data        <- res$input_data
      incProgress(1.0, detail = "Complete!")
    })

    nav_select("sight_main_tabs", selected = "tab_explanation", session = session)
    showNotification("SIGHT screening complete!", type = "message")
  })

  # Reset form
  observeEvent(input$reset_btn, {
    updateTextInput(session, "patient_id", value = "")
    updateNumericInput(session, "total_hgb", value = 14.0)
    updateNumericInput(session, "hematocrit", value = 42.0)
    updateNumericInput(session, "rbc_count", value = 4.8)
    updateNumericInput(session, "mcv", value = 90.0)
    updateNumericInput(session, "mchc", value = 34.0)
    updateNumericInput(session, "rdw", value = 13.0)
    updateSelectInput(session, "threshold_preset", selected = "0.35")

    values$prediction_result <- NULL
    values$shap_data <- NULL
    values$input_data <- NULL
    showNotification("Form reset", type = "message")
  })

  # ---------------- BATCH ANALYSIS HANDLERS ----------------

  output$download_batch_template <- downloadHandler(
    filename = function() {
      "SIGHT_Batch_Template.csv"
    },
    content = function(file) {
      readr::write_csv(get_sight_batch_template(), file)
    }
  )

  execute_sight_batch <- function(raw_df) {
    th_val <- as.numeric(input$batch_threshold_preset %||% "0.35")
    withProgress(message = "Running SIGHT Batch Screening + TreeSHAP...", value = 0, {
      res <- run_sight_batch_analysis(
        raw_df = raw_df,
        model = final_model,
        threshold = th_val,
        progress_callback = function(val, det) {
          setProgress(value = val, detail = det)
        }
      )
      values$batch_results_df <- res$results_df
      values$batch_summary    <- res$summary_stats
      values$batch_error      <- res$error
    })

    if (!is.null(values$batch_error)) {
      showNotification(values$batch_error, type = "error", duration = 8)
    } else {
      nav_select("sight_main_tabs", selected = "tab_batch", session = session)
      showNotification(
        sprintf("Successfully screened %d patients with full SHAP columns!", nrow(values$batch_results_df)),
        type = "message",
        duration = 5
      )
    }
  }

  observeEvent(input$run_batch_btn, {
    req(input$batch_file)
    ext <- tolower(tools::file_ext(input$batch_file$name))
    raw_df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
          stop("Package 'readxl' is required to read Excel files. Please upload a .csv file.")
        }
        readxl::read_excel(input$batch_file$datapath)
      } else {
        readr::read_csv(input$batch_file$datapath, show_col_types = FALSE)
      }
    }, error = function(e) {
      showNotification(paste("Failed to read file:", e$message), type = "error", duration = 6)
      return(NULL)
    })
    req(raw_df)
    execute_sight_batch(raw_df)
  })

  observeEvent(input$load_demo_batch_btn, {
    demo_df <- get_sight_batch_template()
    execute_sight_batch(demo_df)
  })

  output$download_batch_results <- downloadHandler(
    filename = function() {
      paste0("SIGHT_Batch_Predictions_with_SHAP_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$batch_results_df)
      readr::write_csv(values$batch_results_df, file)
    }
  )

  output$download_batch_results_inline <- downloadHandler(
    filename = function() {
      paste0("SIGHT_Batch_Predictions_with_SHAP_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$batch_results_df)
      readr::write_csv(values$batch_results_df, file)
    }
  )

  # ---------------- UI RENDER OUTPUTS ----------------

  output$left_panel_content_ui <- renderUI({
    if (identical(input$analysis_mode, "batch")) {
      render_sight_batch_left_panel()
    } else {
      render_sight_prediction_results(values$prediction_result, values$input_data)
    }
  })

  render_sight_batch_left_panel <- function() {
    format_spec_df <- tibble::tibble(
      `Column Name` = c("total_hgb *", "hematocrit *", "rbc_count *", "mcv *", "mchc *", "rdw *", "patient_id"),
      `Requirement` = c("Required", "Required", "Required", "Required", "Required", "Required", "Optional"),
      `Accepted Aliases` = c(
        "HGB, Hemoglobin, Hb", "HCT, Hematocrit-Blood", "RBC, Red Blood Cell Count",
        "MCV, Mean Cell Volume", "MCHC", "RDW, RDW-CV", "patient_mrn, MRN, ID, sample_id"
      ),
      `Units & Typical Range` = c(
        "g/dL (3.0–20.0)", "% (15.0–60.0)", "M/μL or ×10¹²/L (1.5–8.0)",
        "fL (50.0–120.0)", "g/dL (25.0–40.0)", "% (10.0–35.0)", "Alphanumeric identifier"
      )
    )

    tagList(
      if (!is.null(values$batch_error)) {
        div(class = "alert alert-danger", icon("triangle-exclamation"), " ", strong("Batch Error: "), values$batch_error)
      },

      if (!is.null(values$batch_summary)) {
        s <- values$batch_summary
        div(
          class = "mb-4",
          h5(class = "fw-bold mb-3", icon("chart-pie"), " Cohort Batch Screening Summary"),
          div(
            class = "row g-2",
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", s$total_rows),
              div(class = "batch-kpi-label", "Total Patients Screened")
            )),
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val text-danger", paste0(s$scd_positive, " (", s$positive_rate, "%)")),
              div(class = "batch-kpi-label", "SCD Positive (High Risk)")
            )),
            div(class = "col-md-4", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val text-success", s$scd_negative),
              div(class = "batch-kpi-label", "SCD Negative (Low Risk)")
            )),
            div(class = "col-md-6", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", paste0(s$mean_scd_prob, "%")),
              div(class = "batch-kpi-label", "Mean Predicted SCD Probability")
            )),
            div(class = "col-md-6", div(class = "batch-kpi-box",
              div(class = "batch-kpi-val", paste0(s$threshold_used * 100, "%")),
              div(class = "batch-kpi-label", "Screening Cutoff Applied")
            ))
          ),
          div(
            class = "card mt-2",
            div(class = "card-header py-2", strong("Cohort Mean |SHAP| CBC Feature Importance")),
            div(class = "card-body p-2", plotOutput("batch_shap_summary_plot", height = "250px"))
          )
        )
      } else {
        div(
          class = "alert alert-info mb-3",
          icon("circle-info"),
          strong(" Ready for Batch Screening: "),
          "Upload a CSV or Excel file in the sidebar (or click 'Load Demo Batch') to generate SCD predictions and per-feature TreeSHAP columns."
        )
      },

      div(
        class = "card",
        div(class = "card-header bg-light", strong(icon("list-check"), " Required & Optional Input Spreadsheet Columns")),
        div(
          class = "card-body p-2",
          p(class = "small text-muted mb-2",
            "Column names are matched case-insensitively and common clinical analyzer aliases (e.g., 'HGB', 'HCT', 'RBC') are mapped automatically."),
          DT::renderDataTable({
            datatable(
              format_spec_df,
              options = list(dom = "t", paging = FALSE, ordering = FALSE),
              rownames = FALSE
            )
          })
        )
      )
    )
  }

  output$batch_results_right_ui <- renderUI({
    if (is.null(values$batch_results_df)) {
      return(div(
        class = "alert alert-info",
        icon("table"),
        h4("Batch Screening Results & SHAP Output"),
        p("Switch to 'Batch Spreadsheet Analysis' in the sidebar and upload a CSV/Excel file (or click 'Load Demo Batch') to view row-by-row SCD probabilities and exact TreeSHAP contributions for all 6 CBC parameters.")
      ))
    }

    div(
      div(
        class = "d-flex justify-content-between align-items-center mb-3",
        div(
          h4(class = "mb-0", "Batch Screening Predictions & SHAP Transparency Table"),
          span(class = "text-muted small",
               sprintf("Showing %d rows and %d total columns (including all 6 CBC TreeSHAP columns)",
                       nrow(values$batch_results_df), ncol(values$batch_results_df)))
        ),
        downloadButton("download_batch_results_inline", "Download Full CSV (with SHAP)",
                       class = "btn-success btn-sm")
      ),
      DT::dataTableOutput("batch_results_table")
    )
  })

  output$batch_results_table <- DT::renderDataTable({
    req(values$batch_results_df)
    datatable(
      values$batch_results_df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        autoWidth = TRUE
      ),
      rownames = FALSE
    )
  })

  output$batch_shap_summary_plot <- renderPlot({
    req(values$batch_results_df)
    df_res <- values$batch_results_df
    shap_cols <- c(
      "shap_total_hgb"  = "Total Hemoglobin (g/dL)",
      "shap_hematocrit" = "Hematocrit (%)",
      "shap_rbc_count"  = "RBC Count (M/μL)",
      "shap_mcv"        = "MCV (fL)",
      "shap_mchc"       = "MCHC (g/dL)",
      "shap_rdw"        = "RDW (%)"
    )
    present_cols <- intersect(names(shap_cols), names(df_res))
    req(length(present_cols) > 0)

    imp_tbl <- tibble::tibble(
      Feature = unname(shap_cols[present_cols]),
      MeanAbsSHAP = sapply(present_cols, function(col) mean(abs(as.numeric(df_res[[col]])), na.rm = TRUE))
    ) %>%
      dplyr::arrange(MeanAbsSHAP) %>%
      dplyr::mutate(Feature = factor(Feature, levels = Feature))

    ggplot(imp_tbl, aes(x = Feature, y = MeanAbsSHAP)) +
      geom_col(fill = "#593196", width = 0.65) +
      coord_flip() +
      labs(
        x = NULL,
        y = "Mean |SHAP| Contribution (log-odds)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(face = "bold"),
        panel.grid.major.y = element_blank()
      )
  })

  output$shap_explanation <- renderUI({
    render_sight_shap_explanation(values$shap_data, values$prediction_result)
  })

  output$shap_waterfall <- renderPlot({
    req(values$shap_data)
    values$shap_data$waterfall_plot
  }, res = 96)

  output$shap_contributions_table <- DT::renderDataTable({
    req(values$shap_data)
    contrib_display <- values$shap_data$feature_contributions %>%
      select(
        `CBC Feature` = feature_label_clean,
        `Patient Value` = feature_value,
        `SHAP Value (log-odds)` = shap_value,
        Direction = direction,
        `Clinical Interpretation` = clinical_interpretation
      ) %>%
      mutate(
        `Patient Value` = round(`Patient Value`, 2),
        `SHAP Value (log-odds)` = round(`SHAP Value (log-odds)`, 4)
      )

    datatable(
      contrib_display,
      options = list(
        pageLength = 10,
        dom = "t",
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })

  output$feature_importance_plot <- renderPlotly({
    if (!is.null(values$shap_data) && !is.null(values$shap_data$importance_plot)) {
      ggplotly(values$shap_data$importance_plot, tooltip = c("x", "y"))
    } else {
      plotly_empty() %>%
        layout(
          title = "CBC Feature Importance will appear here after prediction",
          annotations = list(
            text = "Click 'Predict SCD Status' to view interactive SHAP feature importance",
            showarrow = FALSE,
            font = list(size = 15)
          )
        )
    }
  })

  output$clinical_guidelines <- renderUI({
    render_sight_clinical_guidelines()
  })

  output$about_content <- renderUI({
    render_sight_about_content()
  })

  output$download_report <- downloadHandler(
    filename = function() {
      paste0("SIGHT_screening_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      generate_sight_report(values$prediction_result, values$shap_data, values$input_data, file)
    }
  )
}

# ==============================================================================
# 7. RUN APPLICATION
# ==============================================================================

shinyApp(ui = ui, server = server)
