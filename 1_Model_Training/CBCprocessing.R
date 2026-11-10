# SIGHT Data Preparation and Exploratory Data Analysis (EDA)
#
# This script performs the following actions:
# 1. Loads the original 'predictCBC.csv' file.
# 2. Separates the NHANES (controls) and Grady (cases) data.
# 3. Loads and reformats the 'CBC_Ghana.csv' data.
# 4. Loads and reformats the 'CBC_Senegal.csv' data.
# 5. Performs targeted outlier removal on MCHC (hardcoded), RBC count (2 SD rule), and RDW (3 SD rule).
# 6. Saves cleaned, separated datasets for each source.
# 7. Combines all cleaned data into a single, unified dataset ready for model training.
# 8. Generates and saves comparative density plots of the cleaned data.

# === 1. Load Libraries ===
# Ensure you have 'tidyverse' installed: install.packages("tidyverse")
library(tidyverse)
library(patchwork) # For combining plots

cat("--- Starting Data Processing and EDA ---\n\n")

# === 2. Process Original Dataset (predictCBC.csv) ===
tryCatch({
  original_df <- read_csv("predictCBC.csv")
  
  # 2a. Separate and process NHANES data
  nhanes_data_raw <- original_df %>%
    filter(genotype == "AA") %>%
    mutate(
      patient_mrn = as.character(patient_mrn),
      across(c(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw), as.numeric)
    ) %>%
    drop_na(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw)
  
  # 2b. Separate and process Grady data
  grady_data_raw <- original_df %>%
    filter(genotype != "AA") %>%
    mutate(
      patient_mrn = as.character(patient_mrn),
      across(c(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw), as.numeric)
    ) %>%
    drop_na(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw)
  
}, error = function(e) {
  stop("Failed to load or process 'predictCBC.csv'.")
})

# === 3. Process Ghanaian Dataset (CBC_Ghana.csv) ===
tryCatch({
  ghana_raw_df <- read_csv("CBC_Ghana.csv")
  
  ghana_formatted_raw <- ghana_raw_df %>%
    rename(
      total_hgb = HGB, rbc_count = RBC, hematocrit = HCT,
      mcv = `MCV-mean cell volume fl`, mchc = MCHC, rdw = RDW, genotype = GENOTYPE
    ) %>%
    mutate(
      across(c(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw), as.numeric),
      scd_status = case_when(
        genotype %in% c("S-BETA", "SS", "SC", "SF") ~ "SCD",
        genotype %in% c("AA", "AS", "AC", "CC") ~ "No_SCD",
        TRUE ~ NA_character_
      ),
      patient_mrn = paste0("gh_", row_number())
    ) %>%
    select(patient_mrn, genotype, scd_status, hematocrit, mchc, mcv, rbc_count, total_hgb, rdw) %>%
    drop_na(scd_status, hematocrit, mchc, mcv, rbc_count, total_hgb, rdw)
  
}, error = function(e) {
  stop("Failed to load or process 'CBC_Ghana.csv'.")
})

# === 4. Process Senegalese Dataset (CBC_Senegal.csv) ===
tryCatch({
  senegal_raw_df <- read_csv("CBC_Senegal.csv")
  
  senegal_formatted_raw <- senegal_raw_df %>%
    rename(
      patient_mrn = ID, total_hgb = HGB, rbc_count = RBC, hematocrit = HCT,
      mcv = MCV, mchc = MCHC, rdw = RDW
    ) %>%
    mutate(
      patient_mrn = as.character(patient_mrn),
      across(c(hematocrit, mchc, mcv, rbc_count, total_hgb, rdw), as.numeric),
      scd_status = case_when(
        genotype %in% c("SS", "SC","Sbeta+") ~ "SCD",
        genotype %in% c("AA", "AS", "AC") ~ "No_SCD",
        TRUE ~ NA_character_
      )
    ) %>%
    select(patient_mrn, genotype, scd_status, hematocrit, mchc, mcv, rbc_count, total_hgb, rdw) %>%
    drop_na(scd_status, hematocrit, mchc, mcv, rbc_count, total_hgb, rdw)
  
}, error = function(e) {
  stop("Failed to load or process 'CBC_Senegal.csv'.")
})


# === 5. Outlier Identification and Removal (Targeted Rule) ===
cat("\n--- Outlier Analysis and Removal ---\n")

# Function to report and filter outliers using targeted rules
analyze_and_clean <- function(df, name) {
  cat(sprintf("\nAnalyzing dataset: %s (Initial rows: %d)\n", name, nrow(df)))
  
  # Rule 1: Hardcoded physiological range for MCHC
  mchc_lower_hard <- 25
  mchc_upper_hard <- 40
  
  # Rule 2: 2 Standard Deviation rule for RBC Count
  rbc_lower_stat <- mean(df$rbc_count, na.rm = TRUE) - 2 * sd(df$rbc_count, na.rm = TRUE)
  rbc_upper_stat <- mean(df$rbc_count, na.rm = TRUE) + 2 * sd(df$rbc_count, na.rm = TRUE)
  
  # Rule 3: 3 Standard Deviation rule for RDW
  rdw_lower_stat <- mean(df$rdw, na.rm = TRUE) - 3 * sd(df$rdw, na.rm = TRUE)
  rdw_upper_stat <- mean(df$rdw, na.rm = TRUE) + 3 * sd(df$rdw, na.rm = TRUE)
  
  cleaned_df <- df %>%
    filter(
      mchc >= mchc_lower_hard & mchc <= mchc_upper_hard,
      rbc_count >= rbc_lower_stat & rbc_count <= rbc_upper_stat,
      rdw >= rdw_lower_stat & rdw <= rdw_upper_stat
    )
  
  outliers_removed <- nrow(df) - nrow(cleaned_df)
  
  if (outliers_removed > 0) {
    cat(sprintf("Found and removed %d outliers from %s (Targeted rules for MCHC, RBC, and RDW).\n", outliers_removed, name))
  } else {
    cat(sprintf("No outliers found in %s based on targeted rules.\n", name))
  }
  
  return(cleaned_df)
}

# Apply cleaning to all datasets
nhanes_data <- analyze_and_clean(nhanes_data_raw, "NHANES")
grady_data <- analyze_and_clean(grady_data_raw, "Grady")
ghana_formatted <- analyze_and_clean(ghana_formatted_raw, "Ghana")
senegal_formatted <- analyze_and_clean(senegal_formatted_raw, "Senegal")

# === 6. Save Individual Cleaned Data Files ---
write_csv(nhanes_data, "nhanes_controls_cleaned.csv")
cat(sprintf("\n✅ Saved %d cleaned NHANES control samples.\n", nrow(nhanes_data)))
write_csv(grady_data, "grady_scd_cases_cleaned.csv")
cat(sprintf("✅ Saved %d cleaned Grady SCD case samples.\n", nrow(grady_data)))

ghana_controls <- ghana_formatted %>% filter(scd_status == "No_SCD")
write_csv(ghana_controls, "ghana_controls_cleaned.csv")
cat(sprintf("✅ Saved %d cleaned Ghana control samples.\n", nrow(ghana_controls)))
ghana_cases <- ghana_formatted %>% filter(scd_status == "SCD")
write_csv(ghana_cases, "ghana_scd_cases_cleaned.csv")
cat(sprintf("✅ Saved %d cleaned Ghana SCD case samples.\n", nrow(ghana_cases)))

senegal_controls <- senegal_formatted %>% filter(scd_status == "No_SCD")
write_csv(senegal_controls, "senegal_controls_cleaned.csv")
cat(sprintf("✅ Saved %d cleaned Senegal control samples.\n", nrow(senegal_controls)))
senegal_cases <- senegal_formatted %>% filter(scd_status == "SCD")
write_csv(senegal_cases, "senegal_scd_cases_cleaned.csv")
cat(sprintf("✅ Saved %d cleaned Senegal SCD case samples.\n\n", nrow(senegal_cases)))

# === 7. Combine All Cleaned Data for Model Training ===
cat("\n--- Combining All Cleaned Datasets for Training ---\n")

# Add a 'source' column to each dataframe for traceability
combined_training_data <- bind_rows(
  nhanes_data %>% mutate(source = "NHANES"),
  grady_data %>% mutate(source = "Grady"),
  ghana_formatted %>% mutate(source = "Ghana"),
  senegal_formatted %>% mutate(source = "Senegal")
)

# Save the final combined dataset
write_csv(combined_training_data, "combined_cleaned_data_for_training.csv")
cat(sprintf("✅ Saved combined dataset with %d total samples to 'combined_cleaned_data_for_training.csv'.\n\n", nrow(combined_training_data)))

cat("--- Final Dataset Composition ---\n")
print(table(combined_training_data$scd_status, combined_training_data$source))
cat("---------------------------------\n\n")


# === 8. Visualize CBC Parameter Distributions (using CLEANED data) ===
cat("--- Generating CBC Parameter Distribution Plots ---\n")

# Combine dataframes from memory for plotting
combined_plot_df <- bind_rows(
  nhanes_data %>% mutate(source = "NHANES Controls"),
  grady_data %>% mutate(source = "Grady SCD Cases"),
  ghana_controls %>% mutate(source = "Ghana Controls"),
  ghana_cases %>% mutate(source = "Ghana SCD Cases"),
  senegal_controls %>% mutate(source = "Senegal Controls"),
  senegal_cases %>% mutate(source = "Senegal SCD Cases")
) %>%
  mutate(source = factor(source, levels = c(
    "NHANES Controls", "Ghana Controls", "Senegal Controls", 
    "Grady SCD Cases", "Ghana SCD Cases", "Senegal SCD Cases"
  ))) %>%
  select(source, hematocrit, mchc, mcv, rbc_count, total_hgb, rdw)

# Convert to long format for plotting
long_df <- combined_plot_df %>%
  pivot_longer(cols = -source, names_to = "parameter", values_to = "value")

# Create the density plots
distribution_plot <- ggplot(long_df, aes(x = value, fill = source)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  scale_fill_manual(
    values = c(
      "NHANES Controls" = "#1f77b4", "Ghana Controls" = "#aec7e8", "Senegal Controls" = "#729ece",
      "Grady SCD Cases" = "#d62728", "Ghana SCD Cases" = "#ff9896", "Senegal SCD Cases" = "#ff6347"
    )
  ) +
  labs(
    title = "Distribution of CBC Parameters Across Datasets (After Targeted Cleaning)",
    subtitle = "Comparing US, Ghanaian, and Senegalese Populations for SCD Cases and Controls",
    x = "Parameter Value", y = "Density", fill = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# Save the plot
ggsave("cbc_distributions_cleaned_targeted.png", plot = distribution_plot, width = 14, height = 8, dpi = 300)

cat("✅ Successfully saved the cleaned distribution plot to 'cbc_distributions_cleaned_targeted.png'\n")
cat("--- Script Finished ---\n")

print(distribution_plot)

