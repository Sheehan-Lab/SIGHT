================================================================================
SIGHT: Sickle-cell Identification using General Hematology Testing (v2.0)
Standalone Offline Application for Windows
================================================================================

OVERVIEW:
SIGHT is a machine learning-powered clinical decision support tool that screens
for Sickle Cell Disease (SCD) using a single routine Complete Blood Count (CBC).
This standalone executable (SIGHT.exe) runs 100% offline with zero external runtime
dependencies:
  - No R or Python installation required.
  - No internet connection required.
  - Self-contained LightGBM ensemble model (1,583 decision trees) + exact TreeSHAP
    explainability engine embedded directly inside the binary.

--------------------------------------------------------------------------------
HOW TO RUN:
--------------------------------------------------------------------------------
1. Interactive Graphical User Interface (GUI):
   - Simply double-click "SIGHT.exe".
   - A native desktop window will launch along with the clean, modern interactive
     SIGHT user interface.
   - You can choose between:
       * "Single Patient Screening": Enter 6 CBC parameters to get immediate
         SCD probability, screening call, clinical warnings, and an interactive
         TreeSHAP waterfall plot explaining the decision. Download an HTML report.
       * "Batch Spreadsheet Analysis": Upload any CSV or Excel (.xlsx) file, or
         click "Load Demo Batch", to screen an entire cohort and download a
         comprehensive 25-column results spreadsheet including full TreeSHAP values.

2. Direct Batch File Processing (GUI Dialog):
   - Click "Quick Direct Batch Analysis" in the SIGHT control window.
   - Select your input spreadsheet (CSV or Excel) and choose where to save the
     resulting CSV file.

3. Command-Line / Automated Batch Processing:
   - Run SIGHT headlessly from PowerShell or Command Prompt:
       SIGHT.exe --batch input_patients.csv output_results_with_shap.csv [optional_cutoff]
     Example:
       SIGHT.exe --batch cohort.csv results.csv 0.35
   - To generate a blank template CSV from the command line:
       SIGHT.exe --template template.csv

--------------------------------------------------------------------------------
SPREADSHEET FORMAT SPECIFICATIONS:
--------------------------------------------------------------------------------
Required Core CBC Columns (case-insensitive headers):
  1. total_hgb   - Total Hemoglobin (g/dL)
  2. mcv         - Mean Corpuscular Volume (fL)
  3. mchc        - Mean Corpuscular Hemoglobin Concentration (g/dL)
  4. rdw         - Red Cell Distribution Width (%)
  5. hematocrit  - Hematocrit (%)
  6. rbc_count   - Red Blood Cell Count (M/uL or x10^12/L)

Optional Column:
  - patient_id   - Patient identifier or sample code

Output Spreadsheet Columns (25 Total Columns):
  All original columns plus:
  - predicted_status (SCD / No_SCD)
  - screening_call (SCD POSITIVE / SCD NEGATIVE)
  - scd_probability
  - likelihood_of_SCD_percent
  - no_scd_probability
  - threshold_used
  - recommended_action
  - shap_baseline_log_odds
  - shap_prediction_log_odds
  - shap_total_hgb
  - shap_mcv
  - shap_mchc
  - shap_rdw
  - shap_hematocrit
  - shap_rbc_count
  - top_positive_shap_feature
  - top_negative_shap_feature
  - shap_clinical_summary

--------------------------------------------------------------------------------
CLINICAL NOTICE:
--------------------------------------------------------------------------------
SIGHT is intended for research, educational, and frontline clinical decision support.
A positive screening result warrants Stage-2 confirmatory hemoglobin electrophoresis
or HPLC (which can be evaluated using SHARP for definitive 5-class genotyping).
================================================================================
