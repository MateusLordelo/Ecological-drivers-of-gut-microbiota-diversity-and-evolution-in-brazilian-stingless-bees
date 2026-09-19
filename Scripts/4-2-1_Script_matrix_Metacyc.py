import pandas as pd
import numpy as np
from pathlib import Path

# ------------------------------------------------------------------
# USER INPUTS - Only change Base_folder_path
# ------------------------------------------------------------------
Base_folder_path = r"C:\Users\mateu\OneDrive\Area_de_Trabalho\shared_lis\29-spp"

# Convert to a Path object
base = Path(Base_folder_path)

# Define all other paths relative to the base
input_file    = base / "Functional_prediction" / "Metacyc" / "Beta_diversity" / "ANCOM_Full_Total_pairwise.csv"
output_file   = base / "Functional_prediction" / "Metacyc" / "Metacyc_matrix_codes.csv"
metadata_path = base / "Functional_prediction" / "picrust2_out_my_data" / "pathways_out" / "consolidated_pathway_metadata.tsv"
output_file_2 = base / "ASV_tables" / "For_use" / "Metacyc_matrix.csv"

# Handy reference for your script's logic
matrix_input  = output_file


# ------------------------------------------------------------------
# LOADING DATA
# ------------------------------------------------------------------
df = pd.read_csv(input_file)

# We will store data in a "long" format list first:
# [{'Genus': 'Melipona', 'Function': 'PathwayX', 'LFC': 2.5}, ...]
extracted_data = []

# ------------------------------------------------------------------
# PARSING COLUMNS
# ------------------------------------------------------------------
# Find all columns that start with 'lfc_'
lfc_cols = [c for c in df.columns if c.startswith("lfc_")]

# Define your intercept genus here
intercept_genus = "Cephalotrigona"

for col in lfc_cols:
    # 1. Extract Genus names
    parts = col.split("_")
    
    if len(parts) == 3:
        # Standard pairwise format: lfc_GenusA_GenusB
        genus1 = parts[1].replace("Genus", "")
        genus2 = parts[2].replace("Genus", "")
    elif len(parts) == 2:
        # Intercept format: lfc_GenusA
        # In this case, GenusA is compared against the Intercept
        genus1 = parts[1].replace("Genus", "")
        genus2 = intercept_genus
    else:
        print(f"Skipping unexpected column format: {col}")
        continue
    
    # 2. Identify the corresponding significance column
    robust_col = col.replace("lfc_", "diff_robust_")
    
    if robust_col not in df.columns:
        print(f"Warning: Significance column {robust_col} not found. Skipping {col}.")
        continue

    # 3. Iterate through the rows
    for _, row in df.iterrows():
        function_name = row["taxon"]
        lfc_val = row[col]
        
        # Check significance (handling potential NaN or non-bool values)
        is_significant = str(row[robust_col]).strip().upper() == "TRUE"

        if is_significant and pd.notna(lfc_val):
            # Your LFC logic is correct:
            # Positive LFC = Enriched in genus1
            # Negative LFC = Enriched in genus2 (or the Intercept)
            if lfc_val > 0:
                extracted_data.append({
                    "Genus": genus1,
                    "Function": function_name,
                    "LFC_Abs": abs(lfc_val)
                })
            elif lfc_val < 0:
                extracted_data.append({
                    "Genus": genus2,
                    "Function": function_name,
                    "LFC_Abs": abs(lfc_val)
                })
# ------------------------------------------------------------------
# AGGREGATION & MATRIX GENERATION
# ------------------------------------------------------------------
# Convert the list to a DataFrame
long_df = pd.DataFrame(extracted_data)

if long_df.empty:
    print("No significant enrichments found based on the 'TRUE' criteria.")
else:
    # Calculate the mean LFC per Genus and Function
    # This groups by Genus/Function and averages the absolute LFC values
    matrix_df = long_df.groupby(["Genus", "Function"])["LFC_Abs"].mean().unstack()

    # Fill NaN with 0 (where a genus was never enriched in a specific function)
    matrix_df = matrix_df.fillna(0)

    # ------------------------------------------------------------------
    # EXPORT
    # ------------------------------------------------------------------
    # Following your request: Y (index) is Genus, X (columns) is Functions
    matrix_df.to_csv(output_file)
    print(f"Success! Matrix saved to {output_file}")
    print(f"Dimensions: {matrix_df.shape[0]} Genera x {matrix_df.shape[1]} Functions")




import pandas as pd


# ------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------
# Load the matrix (Genus as index, Pathway codes as columns)
matrix_df = pd.read_csv(matrix_input, index_col=0)

# Load the metadata mapping
meta_df = pd.read_csv(metadata_path, sep='\t')

# ------------------------------------------------------------------
# 2. CREATE MAPPING & RENAME
# ------------------------------------------------------------------
# Create a dictionary: {Pathway_Code: Second_Level}
mapping = dict(zip(meta_df['Pathway'], meta_df['Second_Level']))

# Identify columns that have a mapping; keep original name if not found
# (Handy if some codes aren't in your metadata)
new_columns = [mapping.get(col, col) for col in matrix_df.columns]
matrix_df.columns = new_columns

# ------------------------------------------------------------------
# 3. COLLAPSE REPEATED COLUMNS
# ------------------------------------------------------------------
# Group by column name (now Second_Level) and calculate the mean
# axis=1 tells pandas to group columns rather than rows
final_df = matrix_df.groupby(axis=1, level=0).mean()

# ------------------------------------------------------------------
# 4. EXPORT
# ------------------------------------------------------------------
final_df.to_csv(output_file_2)

print(f"Success! Collapsed matrix saved to: {output_file_2}")
print(f"Original pathways: {len(mapping)} mapped to {final_df.shape[1]} unique Second_Level classes.")