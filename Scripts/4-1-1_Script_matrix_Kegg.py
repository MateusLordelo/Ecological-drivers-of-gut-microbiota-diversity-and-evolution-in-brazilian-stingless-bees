import os
import pandas as pd
from pathlib import Path
from collections import defaultdict

# ------------------------------------------------------------------
# USER INPUTS - Only change Base_folder_path
# ------------------------------------------------------------------
Base_folder_path = r"C:/Users/mateu/OneDrive/Area_de_Trabalho/shared_lis/29-spp"

# Convert to Path object
base = Path(Base_folder_path)

# Define directory where all CSV files are stored
input_dir = base / "Functional_prediction" / "Picrust_downstream" / "Kegg_enrichment"

# Define Output file
output_file = base / "ASV_tables" / "For_use" / "Kegg_matrix.csv"

# List of valid groups (used for validation and filtering)
groups = [
    "Cephalotrigona", "Friesella", "Frieseomelitta", "Lestrimelitta",
    "Melipona", "Nannotrigona", "Oxytrigona", "Partamona",
    "Plebeia", "Scaptotrigona", "Schwarziana", "Tetragona",
    "Tetragonisca", "Trigona"
]

# ------------------------------------------------------------------
# DATA STRUCTURE
# ------------------------------------------------------------------

# We will store data as:
# { group: { pathway: [list of LFC values] } }
# defaultdict allows automatic creation of nested dictionaries/lists
enrichment_dict = defaultdict(lambda: defaultdict(list))

# ------------------------------------------------------------------
# FILE PROCESSING LOOP
# ------------------------------------------------------------------

# Iterate over all files in the directory
for filename in os.listdir(input_dir):

    # Only process relevant CSV files
    if not filename.endswith(".csv"):
        continue

    if not filename.startswith("KEGG_Enrichment_LFC_"):
        continue

    # --------------------------------------------------------------
    # PARSE GROUP NAMES FROM FILENAME
    # Expected format:
    # KEGG_Enrichment_LFC_[group1]_vs_[group2]_FullDataset.csv
    # --------------------------------------------------------------
    try:
        parts = filename.replace(".csv", "").split("_")

        # Extract group names based on known pattern
        # ["KEGG", "Enrichment", "LFC", group1, "vs", group2, "FullDataset"]
        group1 = parts[3]
        group2 = parts[5]

    except Exception as e:
        print(f"Skipping file (naming issue): {filename}")
        continue

    # Validate group names
    if group1 not in groups or group2 not in groups:
        print(f"Skipping file (unknown groups): {filename}")
        continue

    # --------------------------------------------------------------
    # READ CSV FILE
    # --------------------------------------------------------------
    file_path = os.path.join(input_dir, filename)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        continue

    # Check required columns exist
    if "mean_LFC" not in df.columns or "PathwayName" not in df.columns:
        print(f"Skipping file (missing columns): {filename}")
        continue

    # --------------------------------------------------------------
    # PROCESS EACH ROW (PATHWAY)
    # --------------------------------------------------------------
    for _, row in df.iterrows():

        pathway = row["PathwayName"]
        lfc = row["mean_LFC"]

        # Skip missing values
        if pd.isna(lfc) or pd.isna(pathway):
            continue

        # ----------------------------------------------------------
        # DETERMINE WHICH GROUP THE PATHWAY IS ENRICHED IN
        # ----------------------------------------------------------
        if lfc > 0:
            # Enriched in group1
            enrichment_dict[group1][pathway].append(lfc)

        elif lfc < 0:
            # Enriched in group2
            enrichment_dict[group2][pathway].append(abs(lfc))
            # Optional: use abs(lfc) so enrichment is always positive magnitude

        # If lfc == 0, ignore (no enrichment)

# ------------------------------------------------------------------
# BUILD FINAL MATRIX
# ------------------------------------------------------------------

# Collect all unique pathways across all groups
all_pathways = set()

for group in enrichment_dict:
    all_pathways.update(enrichment_dict[group].keys())

# Convert to sorted list for consistent ordering
all_pathways = sorted(all_pathways)

# Initialize result DataFrame
result_df = pd.DataFrame(0, index=groups, columns=all_pathways)

# ------------------------------------------------------------------
# FILL MATRIX WITH MEAN VALUES
# ------------------------------------------------------------------

for group in enrichment_dict:
    for pathway in enrichment_dict[group]:

        values = enrichment_dict[group][pathway]

        if len(values) > 0:
            mean_value = sum(values) / len(values)
            result_df.loc[group, pathway] = mean_value

# ------------------------------------------------------------------
# EXPORT RESULT
# ------------------------------------------------------------------

result_df.to_csv(output_file)

print("Enrichment matrix successfully saved to:", output_file)