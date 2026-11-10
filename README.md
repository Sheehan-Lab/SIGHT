# SIGHT: Sickle-cell Identification using General Hematology Testing

SIGHT is a machine learning-based screening tool for Sickle Cell Disease (SCD) that uses standard parameters from a routine Complete Blood Count (CBC). Developed by the Sheehan Lab, SIGHT provides early risk stratification in settings where hemoglobin electrophoresis or high-performance liquid chromatography (HPLC) is unavailable, delayed, or cost-prohibitive.

The underlying model is a LightGBM ensemble trained on routine CBC indices, paired with TreeSHAP feature attributions for transparent, case-by-case explainability.

## Repository Structure

- `1_Model_Training/`
  Code and workflows for data preprocessing, feature engineering, and model training in R (`sight_v6.R`, `CBCprocessing.R`), along with trained LightGBM model objects (`scd_final_model.rds`).
- `2_Shiny_App/`
  Interactive R Shiny web application (`app.R`) for deployment to local environments or hosting platforms such as shinyapps.io. Includes single-sample screening and batch CSV upload with SHAP visualization.
- `3_Standalone_Executable/`
  Self-contained Windows executable (`SIGHT.exe`), batch template (`SIGHT_Batch_Template.csv`), and instructions for offline use without requiring R, Python, or an internet connection.

## Quick Start

### 1. Standalone Windows Executable
If you are on Windows and want to run SIGHT locally without installing R:
1. Navigate to `3_Standalone_Executable/` and double-click `SIGHT.exe`.
2. Enter values manually in the Single Patient tab, or upload a formatted spreadsheet in the Batch tab.
3. For command-line batch screening:
   ```cmd
   SIGHT.exe --batch input_patients.csv output_results.csv
   ```

### 2. R Shiny Application
To run the web interface locally in R:
```R
install.packages(c("shiny", "bslib", "DT", "lightgbm", "shapr", "ggplot2", "plotly", "readr", "readxl"))
shiny::runApp("2_Shiny_App")
```

### 3. Model Training Pipeline
The training pipeline requires R and the packages listed in `1_Model_Training/sight_v6.R`. Training scripts accept de-identified CBC datasets structured according to the standard feature schema.

## Input Features

SIGHT accepts the following CBC indices:
- `total_hgb`: Total Hemoglobin (g/dL)
- `mcv`: Mean Corpuscular Volume (fL)
- `mchc`: Mean Corpuscular Hemoglobin Concentration (g/dL)
- `rdw`: Red Cell Distribution Width (%)
- `hematocrit`: Hematocrit (%)
- `rbc_count`: Red Blood Cell Count (x10^6/uL)
- `patient_id`: Optional identifier string

A formatted sample file is available at `3_Standalone_Executable/SIGHT_Batch_Template.csv`.

## Output and Interpretation

- **SCD Probability**: Predicted likelihood of sickle cell disease (0.0 to 1.0).
- **Screening Call**: Flagged as Positive (high risk) or Negative based on calibrated clinical decision thresholds.
- **TreeSHAP Attributions**: Feature-level contribution values indicating how each individual CBC parameter influenced the final risk score.

## License

This project is released under a dual noncommercial licensing structure:
- **Software Code**: The source code, application scripts, user interface assets, and standalone executables are licensed under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/).
- **Models and Data**: The trained model weights, parameters, serialized workflow files (`.rds`), and documentation are licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)](https://creativecommons.org/licenses/by-nc-sa/4.0/) license.

Commercial use of either the software code or the model artifacts is strictly prohibited without explicit written permission from the Sheehan Lab.

## Citation

If you use SIGHT in your research or clinical workflows, please cite the corresponding publication from the Sheehan Lab.
