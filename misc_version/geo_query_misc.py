#!/usr/bin/env Rscript

############################################################

# Install GEOparse
pip install GEOparse

############################################################

import GEOparse
import os
import shutil
import pandas as pd
import glob
import gzip


#pwd='/scratch16/abattle4/surya/datasets/for_diptavo/downloads/final'

# Download and load the dataset GSE188486
gse = GEOparse.get_GEO(geo="GSE188486")
gse = GEOparse.get_GEO(geo="GSE188486", destdir="./")

# Print the metadata
print("Metadata:")
print(gse.metadata)

# Downloading all supplementary files associated with GSE188486
for supp_file in gse.supplementary_files:
    GEOparse.utils.download_from_url(supp_file, filename=supp_file.split('/')[-1])

metadata_df = gse.phenotype_data
meta_df.to_csv("freedman_etal_metadata.tsv", sep="\t", header=True, index=True)

# Selected GSE dataframe
gse_df_select = gse.phenotype_data[["title", "geo_accession", 'source_name_ch1', 'description', 'supplementary_file_2']]
gse_df_select.to_csv("Freedman_etal_metadata_selectedCols.tsv", sep="\t", header=True, index=True)

# Example of supplemental download file
gse_supp = gse.phenotype_data['supplementary_file_2'][0]

# Filter Histone H3K PTMs
filt_df = gse_df[gse_df["title"].str.contains("H3K")]

# Filter and print unique items from a DataFrame (if applicable)
# Assuming filt_df is your DataFrame and "source_name_ch1" is the column of interest
unique_values = filt_df["source_name_ch1"].unique()
print(unique_values)

# Create a directory for renamed files
os.makedirs('renamed_files', exist_ok=True)

# Iterate over the rows of the dataframe
for index, row in filt_df.iterrows():
    # Extract geo_accession, title, and supplementary file URL
    geo_accession = row['geo_accession']
    title = row['title']
    source_name = row['source_name_ch1'].replace(' ', '_')  # Replace spaces with dashes
    supp_file_url = row['supplementary_file_2']
    
    # Extract the filename from the URL
    filename = supp_file_url.split('/')[-1] if pd.notna(supp_file_url) else None
    print(f"Processing file: {filename}")
    
    # Check if the file exists in the current directory
    if filename and os.path.exists(filename):
        # Create a new filename with source name ID
        new_filename = f"{source_name}-{filename}"
        #new_filename = f"{geo_accession}_{title}_{source_name}_{filename}"
        print(f"Renamed file {filename} to {new_filename}\n")
        
        # Copy and rename the file to the new directory
        shutil.copy(filename, os.path.join('renamed_files', new_filename))

