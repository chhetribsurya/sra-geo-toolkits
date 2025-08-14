# Genomics Data (SRA/GEO) Processing Pipeline

A suite of scripts for downloading and processing genomics data from NCBI SRA and GEO databases. This pipeline provides robust tools for acquiring sequencing data (SRA) and gene expression datasets (GEO) with automated processing, validation, and organization.

## Table of Contents

- [Overview](#overview)
- [Scripts Description](#scripts-description)
- [Installation](#installation)
- [SRA Download Script](#sra-download-script)
- [GEO Data Processing (R)](#geo-data-processing-r)
- [GEO Data Processing (Python)](#geo-data-processing-python)
- [Configuration](#configuration)
- [Troubleshooting](#troubleshooting)
- [Best Practices](#best-practices)

## Overview

This pipeline consists of three main components designed to handle different aspects of genomics data acquisition and processing:

1. **SRA Download Script** (`sra_script.sh`) - Shell-based pipeline for downloading and converting NCBI SRA sequencing data
2. **GEO R Script** (`geo_query.R`) - R-based tool for downloading GEO datasets with Bioconductor integration
3. **GEO Python Script** (`geo_processor.py`) - Python-based comprehensive GEO data management system

Each script is designed to work independently or as part of a larger analysis workflow, providing flexibility for different research needs and computational environments.

## Scripts Description

### SRA Download Script
Automates the complete workflow from SRA experiment identification through FASTQ file preparation. Handles batch processing, file validation, and format conversion with multi-threading support.

### GEO R Script  
Integrates with the Bioconductor ecosystem to download GEO datasets, extract expression matrices, and prepare data for statistical analysis within R environments.

### GEO Python Script
Provides enterprise-level data management for GEO datasets with advanced file organization, metadata processing, and batch workflow capabilities.

## Installation

### System Requirements

**For SRA Script:**
- Unix-like operating system (Linux, macOS)
- Conda package manager (recommended)
- At least 4GB RAM for typical datasets
- Sufficient storage space (SRA files can be large)

**For R Script:**
- R version 4.0 or higher
- Internet connection for package installation
- 8GB RAM recommended for large datasets

**For Python Script:**
- Python 3.7 or higher
- pip package manager
- 8GB RAM recommended for large datasets

### Dependencies Installation

**SRA Script Dependencies:**
```bash
# Using conda (recommended)
conda install -c bioconda entrez-direct sra-tools parallel-fastq-dump pigz

# Or install individually
conda install entrez-direct
conda install sra-tools
conda install -c bioconda parallel-fastq-dump
conda install pigz
```

**R Script Dependencies:**
```r
# Install within R console
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GEOquery")

# Additional recommended packages
install.packages(c("dplyr", "readr", "stringr"))
```

**Python Script Dependencies:**
```bash
pip install GEOparse pandas numpy

# Optional but recommended
pip install jupyter matplotlib seaborn
```

### Script Preparation

1. Clone or download the scripts to your working directory
2. Make the shell script executable:
```bash
chmod +x sra_script.sh
```
3. Verify installations by running help commands:
```bash
./sra_script.sh --help
```

## SRA Download Script

### Basic Usage

The SRA script provides a command-line interface for downloading and processing sequencing data from NCBI's Sequence Read Archive.

**Install Dependencies:**
```bash
./sra_script.sh install
```

**Download Single Dataset:**
```bash
# Complete workflow from SRX to FASTQ
./sra_script.sh download SRX7805363

# Download specific SRR file
./sra_script.sh download SRR11184871
```

**Batch Processing:**
```bash
# Create a file containing SRR IDs (one per line)
echo -e "SRR11184871\nSRR11184872\nSRR11184873" > sra_list.txt

# Process the entire list
./sra_script.sh download-list sra_list.txt
```

### Advanced Options

**Custom Configuration:**
```bash
# Specify threads and output directory
./sra_script.sh -t 16 -o /path/to/output download SRX7805363

# Convert existing SRA file to FASTQ
./sra_script.sh -t 8 convert SRR11184871
```

**Validation and Verification:**
```bash
# Check file integrity
./sra_script.sh validate SRR11184871

# Determine library layout (single vs paired-end)
./sra_script.sh check-layout SRR11184871
```

### Output Structure

```
output_directory/
├── SRR11184871.sra              # Downloaded SRA file
├── fastq_output/
│   ├── SRR11184871_1.fastq.gz   # Forward reads (paired-end)
│   ├── SRR11184871_2.fastq.gz   # Reverse reads (paired-end)
│   └── SRR11184871.fastq.gz     # Single reads (single-end)
└── logs/
    └── download.log             # Processing log
```

### Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-t, --threads` | Number of processing threads | 10 |
| `-o, --output` | Output directory | ./fastq_output |
| `-h, --help` | Show help message | - |

## GEO Data Processing (R)

### Basic Usage

The R script provides functions for downloading and processing GEO datasets within the R environment.

**Load and Initialize:**
```r
source("geo_query.R")

# This will automatically install required packages
```

**Simple Analysis:**
```r
# Download and analyze single dataset
result <- main_geo_analysis("GSE188486", "./geo_data")

# With sample filtering
result <- main_geo_analysis("GSE188486", "./geo_data", "H3K")
```

**Manual Workflow:**
```r
# Step-by-step processing
install_geo_packages()
dataset <- download_geo_dataset("GSE188486", "./manual_analysis")
expr_data <- extract_expression_data(dataset$gse_object, log_transform = TRUE)
metadata <- extract_metadata(dataset$gse_object)
```

### Advanced Features

**Batch Processing:**
```r
# Process multiple datasets
gse_list <- c("GSE188486", "GSE123456", "GSE789012")
batch_results <- process_geo_datasets_batch(gse_list, "./batch_geo_data")
```

**Metadata Filtering:**
```r
# Filter samples by criteria
filtered_metadata <- filter_samples_by_metadata(
    metadata, 
    filter_column = "title", 
    filter_pattern = "H3K27ac",
    case_sensitive = FALSE
)
```

**Expression Data Processing:**
```r
# Extract and transform expression data
expr_data <- extract_expression_data(gse_object, log_transform = TRUE)

# Create sample annotation
annotation <- create_sample_annotation(
    metadata, 
    key_columns = c("title", "geo_accession", "source_name_ch1")
)
```

### Output Structure

```
geo_analysis/
├── GSE188486_metadata.txt       # Sample metadata
├── GSE188486_summary.txt        # Dataset summary
├── GSE188486_annotation.txt     # Sample annotation
└── supplementary_files/         # Downloaded supplementary files
    ├── GSM123_peaks.bed.gz
    └── GSM124_peaks.bed.gz
```

### Function Reference

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `download_geo_dataset()` | Download GEO dataset | gse_id, output_dir, get_supplementary |
| `extract_expression_data()` | Get expression matrix | gse_object, log_transform |
| `filter_samples_by_metadata()` | Filter samples | metadata, filter_column, filter_pattern |
| `generate_dataset_summary()` | Create summary stats | gse_data, output_file |

## GEO Data Processing (Python)

### Basic Usage

The Python script provides a comprehensive command-line interface for GEO data management.

**Complete Analysis Workflow:**
```bash
# Analyze single dataset with full processing
python geo_processor.py analyze GSE188486 --output-dir ./results

# Filter samples during analysis
python geo_processor.py analyze GSE188486 \
    --filter-pattern "H3K27ac" \
    --filter-column "title" \
    --output-dir ./filtered_results
```

**Download Only:**
```bash
# Download dataset without processing
python geo_processor.py download GSE188486 --supplementary

# Download to specific directory
python geo_processor.py download GSE188486 \
    --output-dir ./raw_downloads \
    --supplementary
```

**Batch Processing:**
```bash
# Process multiple datasets
python geo_processor.py batch GSE188486 GSE123456 GSE789012 \
    --output-dir ./batch_results \
    --filter-pattern "ChIP-seq"
```

### Advanced Configuration

**Selective Metadata Extraction:**
```bash
# Extract specific metadata columns
python geo_processor.py analyze GSE188486 \
    --selected-columns title geo_accession source_name_ch1 \
    --output-dir ./selective_analysis
```

**Custom Filtering:**
```bash
# Filter on different columns with regex patterns
python geo_processor.py analyze GSE188486 \
    --filter-pattern ".*kidney.*" \
    --filter-column "source_name_ch1" \
    --output-dir ./kidney_samples
```

### Programming Interface

**Using as Python Module:**
```python
from geo_processor import GEODataProcessor, analyze_geo_dataset

# Initialize processor
processor = GEODataProcessor("./analysis_output")

# Download and process dataset
gse = processor.download_geo_dataset("GSE188486")
metadata = processor.extract_metadata(gse, output_file="metadata.tsv")

# Filter samples
filtered = processor.filter_samples_by_criteria(
    metadata, "title", "H3K27ac", case_sensitive=False
)

# Generate comprehensive summary
summary = processor.generate_dataset_summary(gse, metadata)
```

### Output Structure

```
geo_analysis/
├── GSE188486/
│   ├── GSE188486_metadata.tsv           # Complete metadata
│   ├── GSE188486_filtered_metadata.tsv  # Filtered samples
│   ├── GSE188486_summary.json           # Rich dataset summary
│   ├── GSE188486_sample_annotation.tsv  # Clean annotation
│   ├── renamed_files/                   # Organized supplementary files
│   │   ├── HEK293_Cell_Line-GSM5567890_peaks.bed.gz
│   │   ├── Kidney_Tissue-GSM5567891_peaks.bed.gz
│   │   └── ...
│   └── logs/                           # Processing logs
├── geo_analysis.log                    # Main processing log
└── batch_summary.json                 # Batch processing results
```

### Class Reference

| Class | Purpose | Key Methods |
|-------|---------|-------------|
| `GEODataProcessor` | Main processing class | download_geo_dataset, extract_metadata, filter_samples |

## Configuration

### Environment Variables

You can set environment variables to configure default behavior:

```bash
# SRA Script
export SRA_OUTPUT_DIR="/data/sra_downloads"
export SRA_THREADS=16

# Python Script
export GEO_OUTPUT_DIR="/data/geo_analysis"
export GEO_LOG_LEVEL="DEBUG"
```

### Configuration Files

**Python Script Configuration (config.json):**
```json
{
    "download_supplementary": true,
    "filter_pattern": "H3K",
    "filter_column": "title",
    "selected_columns": ["title", "geo_accession", "source_name_ch1"],
    "output_format": "tsv",
    "compression": true
}
```

Use with:
```bash
python geo_processor.py workflow GSE188486 --config config.json
```

## Troubleshooting

### Common Issues

**SRA Script Issues:**

*Problem: "Command not found" errors*
- Solution: Ensure all dependencies are installed and in PATH
- Check: `which prefetch` and `which fasterq-dump`

*Problem: Download timeouts*
- Solution: Increase timeout settings or try again later
- Use: `prefetch --max-size 50G` for large files

*Problem: Insufficient disk space*
- Solution: SRA files can be very large; ensure adequate storage
- Monitor: Use `df -h` to check available space

**R Script Issues:**

*Problem: Package installation failures*
- Solution: Update R and Bioconductor versions
- Try: `BiocManager::install(version = "3.18")`

*Problem: Memory issues with large datasets*
- Solution: Increase R memory limit
- Use: `memory.limit(size=16000)` on Windows

**Python Script Issues:**

*Problem: GEOparse installation issues*
- Solution: Use conda instead of pip
- Try: `conda install -c bioconda geoparse`

*Problem: Permission denied errors*
- Solution: Check directory permissions
- Use: `chmod 755 output_directory`

### Performance Optimization

**For Large Datasets:**
- Use SSD storage for better I/O performance
- Increase thread counts on multi-core systems
- Process datasets in smaller batches to manage memory

**Memory Management:**
- Monitor memory usage with `htop` or `top`
- For R: Use `gc()` to force garbage collection
- For Python: Process datasets sequentially for large batches

### Logging and Debugging

**Enable Verbose Logging:**
```bash
# SRA Script - check log files in output directory
tail -f sra_download.log

# Python Script - set debug level
python geo_processor.py analyze GSE188486 --log-level DEBUG
```

**Debug Mode:**
```r
# R Script - enable detailed output
options(verbose = TRUE)
```

### Workflow Integration

**Typical Analysis Workflow:**
1. Use SRA script to download sequencing data
2. Use Python GEO script to organize metadata and supplementary files
3. Use R GEO script to query GEO/GSE, organize metadata and supplementary files

### Security and Ethics

1. **Respect data usage policies** of public repositories
2. **Follow institutional guidelines** for data handling
3. **Ensure proper attribution** when using public datasets
4. **Consider data privacy implications** when sharing results
