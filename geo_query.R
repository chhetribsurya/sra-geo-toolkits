#!/usr/bin/env Rscript

# GEO Dataset Download and Analysis Script using GEOquery
# Dependencies: BiocManager, GEOquery

#=============================================================================
# PACKAGE INSTALLATION AND LOADING
#=============================================================================

#' Install required Bioconductor packages
#' @description Installs BiocManager and GEOquery if not already present
install_geo_packages <- function() {
  cat("Installing required packages...\n")
  
  # Install BiocManager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    cat("BiocManager installed successfully.\n")
  }
  
  # Install GEOquery
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery")
    cat("GEOquery installed successfully.\n")
  }
  
  # Load required libraries
  library(GEOquery)
  cat("All packages loaded successfully.\n")
}

#=============================================================================
# GEO DATA DOWNLOAD FUNCTIONS
#=============================================================================

#' Download GEO dataset by GSE accession
#' @param gse_id Character string of GSE accession number (e.g., "GSE188486")
#' @param output_dir Character string specifying output directory (default: current directory)
#' @param get_supplementary Logical indicating whether to download supplementary files (default: TRUE)
#' @return List containing GEO dataset object and metadata
#' @examples
#' dataset <- download_geo_dataset("GSE188486", "./geo_data")
download_geo_dataset <- function(gse_id, output_dir = ".", get_supplementary = TRUE) {
  
  cat(sprintf("Downloading GEO dataset: %s\n", gse_id))
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }
  
  # Set working directory temporarily
  original_wd <- getwd()
  setwd(output_dir)
  
  tryCatch({
    # Download GEO dataset
    gse <- getGEO(gse_id, GSEMatrix = TRUE)
    cat(sprintf("Successfully downloaded dataset: %s\n", gse_id))
    
    # Handle multiple series
    if (length(gse) > 1) {
      cat("Multiple series detected. Using the first series.\n")
      gse_data <- gse[[1]]
    } else {
      gse_data <- gse[[1]]
    }
    
    # Extract metadata
    metadata <- pData(phenoData(gse_data))
    cat(sprintf("Extracted metadata for %d samples\n", nrow(metadata)))
    
    # Download supplementary files if requested
    supp_files <- NULL
    if (get_supplementary && length(gse_data@supplementary_file) > 0) {
      cat("Downloading supplementary files...\n")
      supp_files <- download_supplementary_files(gse_data)
    }
    
    # Reset working directory
    setwd(original_wd)
    
    return(list(
      gse_object = gse_data,
      metadata = metadata,
      supplementary_files = supp_files,
      gse_id = gse_id
    ))
    
  }, error = function(e) {
    setwd(original_wd)
    stop(sprintf("Error downloading dataset %s: %s", gse_id, e$message))
  })
}

#' Download supplementary files for a GEO dataset
#' @param gse_object GEO dataset object
#' @return Character vector of downloaded file names
download_supplementary_files <- function(gse_object) {
  
  if (length(gse_object@supplementary_file) == 0) {
    cat("No supplementary files available.\n")
    return(NULL)
  }
  
  downloaded_files <- character()
  
  for (file_url in gse_object@supplementary_file) {
    tryCatch({
      file_name <- basename(file_url)
      download.file(file_url, destfile = file_name, mode = "wb")
      downloaded_files <- c(downloaded_files, file_name)
      cat(sprintf("Downloaded: %s\n", file_name))
    }, error = function(e) {
      cat(sprintf("Failed to download %s: %s\n", file_url, e$message))
    })
  }
  
  cat(sprintf("Downloaded %d supplementary files\n", length(downloaded_files)))
  return(downloaded_files)
}

#=============================================================================
# DATA EXTRACTION AND PROCESSING FUNCTIONS
#=============================================================================

#' Extract expression data from GEO dataset
#' @param gse_object GEO dataset object
#' @param log_transform Logical indicating whether to log-transform data (default: FALSE)
#' @return Matrix of expression data
extract_expression_data <- function(gse_object, log_transform = FALSE) {
  
  cat("Extracting expression data...\n")
  
  # Extract expression matrix
  expr_data <- exprs(gse_object)
  
  if (log_transform) {
    cat("Applying log2 transformation...\n")
    expr_data <- log2(expr_data + 1)  # Add 1 to avoid log(0)
  }
  
  cat(sprintf("Expression data dimensions: %d genes x %d samples\n", 
              nrow(expr_data), ncol(expr_data)))
  
  return(expr_data)
}

#' Extract and process metadata from GEO dataset
#' @param gse_object GEO dataset object
#' @param selected_columns Character vector of column names to extract (optional)
#' @param output_file Character string for output file name (optional)
#' @return Data frame containing processed metadata
extract_metadata <- function(gse_object, selected_columns = NULL, output_file = NULL) {
  
  cat("Extracting metadata...\n")
  
  # Extract metadata
  metadata <- pData(phenoData(gse_object))
  
  # Select specific columns if specified
  if (!is.null(selected_columns)) {
    available_cols <- intersect(selected_columns, colnames(metadata))
    if (length(available_cols) > 0) {
      metadata <- metadata[, available_cols, drop = FALSE]
      cat(sprintf("Selected %d columns from metadata\n", length(available_cols)))
    } else {
      cat("Warning: None of the specified columns found in metadata\n")
    }
  }
  
  # Save to file if specified
  if (!is.null(output_file)) {
    write.table(metadata, file = output_file, sep = "\t", 
                row.names = TRUE, col.names = TRUE, quote = FALSE)
    cat(sprintf("Metadata saved to: %s\n", output_file))
  }
  
  return(metadata)
}

#' Filter samples based on metadata criteria
#' @param metadata Data frame containing sample metadata
#' @param filter_column Character string specifying column to filter on
#' @param filter_pattern Character string or regex pattern for filtering
#' @param case_sensitive Logical indicating case sensitivity (default: FALSE)
#' @return Data frame with filtered samples
filter_samples_by_metadata <- function(metadata, filter_column, filter_pattern, case_sensitive = FALSE) {
  
  if (!filter_column %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", filter_column))
  }
  
  cat(sprintf("Filtering samples based on '%s' column with pattern: %s\n", 
              filter_column, filter_pattern))
  
  # Apply filter
  if (case_sensitive) {
    filtered_metadata <- metadata[grepl(filter_pattern, metadata[[filter_column]]), ]
  } else {
    filtered_metadata <- metadata[grepl(filter_pattern, metadata[[filter_column]], ignore.case = TRUE), ]
  }
  
  cat(sprintf("Filtered from %d to %d samples\n", 
              nrow(metadata), nrow(filtered_metadata)))
  
  return(filtered_metadata)
}

#=============================================================================
# ANALYSIS AND SUMMARY FUNCTIONS
#=============================================================================

#' Generate summary statistics for GEO dataset
#' @param gse_data List containing GEO dataset information
#' @param output_file Character string for output file name (optional)
#' @return Data frame containing summary statistics
generate_dataset_summary <- function(gse_data, output_file = NULL) {
  
  cat("Generating dataset summary...\n")
  
  gse_object <- gse_data$gse_object
  metadata <- gse_data$metadata
  
  # Basic statistics
  summary_stats <- data.frame(
    Dataset_ID = gse_data$gse_id,
    Total_Samples = nrow(metadata),
    Total_Features = nrow(exprs(gse_object)),
    Platform = annotation(gse_object),
    Title = gse_object@header$title,
    Summary = gse_object@header$summary,
    Submission_Date = gse_object@header$submission_date,
    Last_Update = gse_object@header$last_update_date,
    stringsAsFactors = FALSE
  )
  
  # Sample type distribution if available
  if ("source_name_ch1" %in% colnames(metadata)) {
    sample_types <- table(metadata$source_name_ch1)
    cat("Sample type distribution:\n")
    print(sample_types)
  }
  
  # Save summary if specified
  if (!is.null(output_file)) {
    write.table(summary_stats, file = output_file, sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat(sprintf("Summary saved to: %s\n", output_file))
  }
  
  return(summary_stats)
}

#' Create sample annotation file
#' @param metadata Data frame containing sample metadata
#' @param key_columns Character vector of important columns to include
#' @param output_file Character string for output file name
create_sample_annotation <- function(metadata, key_columns = NULL, output_file = "sample_annotation.txt") {
  
  cat("Creating sample annotation file...\n")
  
  # Default key columns if not specified
  if (is.null(key_columns)) {
    key_columns <- c("title", "geo_accession", "source_name_ch1", 
                     "characteristics_ch1", "description")
  }
  
  # Select available key columns
  available_cols <- intersect(key_columns, colnames(metadata))
  
  if (length(available_cols) > 0) {
    annotation <- metadata[, available_cols, drop = FALSE]
  } else {
    annotation <- metadata
  }
  
  # Save annotation file
  write.table(annotation, file = output_file, sep = "\t", 
              row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  cat(sprintf("Sample annotation saved to: %s\n", output_file))
  return(annotation)
}

#=============================================================================
# BATCH PROCESSING FUNCTIONS
#=============================================================================

#' Process multiple GEO datasets in batch
#' @param gse_ids Character vector of GSE accession numbers
#' @param output_base_dir Character string for base output directory
#' @param get_supplementary Logical indicating whether to download supplementary files
#' @return List containing results for each dataset
process_geo_datasets_batch <- function(gse_ids, output_base_dir = "./geo_batch", get_supplementary = TRUE) {
  
  cat(sprintf("Processing %d GEO datasets in batch...\n", length(gse_ids)))
  
  # Create base output directory
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  results <- list()
  
  for (gse_id in gse_ids) {
    cat(sprintf("\n--- Processing %s ---\n", gse_id))
    
    # Create individual output directory
    dataset_dir <- file.path(output_base_dir, gse_id)
    
    tryCatch({
      # Download dataset
      dataset_result <- download_geo_dataset(gse_id, dataset_dir, get_supplementary)
      
      # Generate summary
      summary_file <- file.path(dataset_dir, paste0(gse_id, "_summary.txt"))
      summary_stats <- generate_dataset_summary(dataset_result, summary_file)
      
      # Create annotation file
      annotation_file <- file.path(dataset_dir, paste0(gse_id, "_annotation.txt"))
      annotation <- create_sample_annotation(dataset_result$metadata, output_file = annotation_file)
      
      # Store results
      results[[gse_id]] <- list(
        dataset = dataset_result,
        summary = summary_stats,
        annotation = annotation,
        output_dir = dataset_dir
      )
      
      cat(sprintf("Successfully processed %s\n", gse_id))
      
    }, error = function(e) {
      cat(sprintf("Error processing %s: %s\n", gse_id, e$message))
      results[[gse_id]] <- list(error = e$message)
    })
  }
  
  cat(sprintf("\nBatch processing completed. Results saved in: %s\n", output_base_dir))
  return(results)
}

#=============================================================================
# MAIN EXECUTION FUNCTION
#=============================================================================

#' Main function to demonstrate usage
#' @param gse_id Character string of GSE accession number
#' @param output_dir Character string for output directory
#' @param filter_pattern Character string for sample filtering (optional)
main_geo_analysis <- function(gse_id = "GSE188486", output_dir = "./geo_analysis", filter_pattern = NULL) {
  
  cat("=== GEO Dataset Analysis Pipeline ===\n")
  
  # Install packages if needed
  install_geo_packages()
  
  # Download dataset
  cat("\n1. Downloading dataset...\n")
  dataset <- download_geo_dataset(gse_id, output_dir)
  
  # Extract data (not required for metadata extraction, but just in case 
  # interested in expression data as well)
  cat("\n2. Extracting expression data (Note:OPTIONAL)...\n")
  expr_data <- extract_expression_data(dataset$gse_object)
  
  # Process metadata
  cat("\n3. Processing metadata...\n")
  metadata_file <- file.path(output_dir, paste0(gse_id, "_metadata.txt"))
  metadata <- extract_metadata(dataset$gse_object, output_file = metadata_file)
  
  # Filter samples if pattern provided
  if (!is.null(filter_pattern)) {
    cat("\n4. Filtering samples...\n")
    filtered_metadata <- filter_samples_by_metadata(metadata, "title", filter_pattern)
    
    # Save filtered metadata
    filtered_file <- file.path(output_dir, paste0(gse_id, "_filtered_metadata.txt"))
    write.table(filtered_metadata, file = filtered_file, sep = "\t", 
                row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  
  # Generate summary
  cat("\n5. Generating summary...\n")
  summary_file <- file.path(output_dir, paste0(gse_id, "_summary.txt"))
  summary_stats <- generate_dataset_summary(dataset, summary_file)
  
  cat("\n=== Analysis completed successfully ===\n")
  cat(sprintf("Results saved in: %s\n", output_dir))
  
  return(list(
    dataset = dataset,
    expression_data = expr_data,
    metadata = metadata,
    summary = summary_stats
  ))
}

#=============================================================================
# USAGE EXAMPLES
#=============================================================================

# Example 1: Download single dataset
# result <- main_geo_analysis("GSE188486", "./geo_data")

# Example 2: Download dataset with filtering
# result <- main_geo_analysis("GSE188486", "./geo_data", "H3K")

# Example 3: Batch processing
# gse_list <- c("GSE188486", "GSE123456", "GSE789012")
# batch_results <- process_geo_datasets_batch(gse_list, "./batch_geo_data")

# Example 4: Manual workflow
# install_geo_packages()
# dataset <- download_geo_dataset("GSE188486", "./manual_analysis")
# expr_data <- extract_expression_data(dataset$gse_object, log_transform = TRUE) #Note: Not reqd. for metadata
# metadata <- extract_metadata(dataset$gse_object, 
#                              selected_columns = c("title", "geo_accession", "source_name_ch1"))
# filtered_metadata <- filter_samples_by_metadata(metadata, "title", "H3K", case_sensitive = FALSE)

cat("GEO Query R Script loaded successfully.\n")
cat("Use main_geo_analysis() to start analysis or call individual functions as needed.\n")
