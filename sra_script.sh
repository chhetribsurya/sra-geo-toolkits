#!/bin/bash

# SRA Download and Processing Script
# Description: Script for downloading and processing SRA files
# Dependencies: entrez-direct, sra-toolkit, pigz, parallel-fastq-dump

set -euo pipefail  # Exit on error, undefined vars, pipe failures

#=============================================================================
# CONFIGURATION
#=============================================================================

# Default settings
DEFAULT_THREADS=10
DEFAULT_OUTPUT_DIR="./fastq_output"
SRA_TOOLS_PATH="sratoolkit.3.0.1-mac64/bin"

#=============================================================================
# UTILITY FUNCTIONS
#=============================================================================

# Print usage information
usage() {
    cat << EOF
Usage: $0 [OPTIONS] COMMAND

COMMANDS:
    install                 Install required dependencies
    download SRX_ID         Download SRA files from SRX ID
    download-list FILE      Download multiple SRA files from list
    convert SRR_ID          Convert SRA to FASTQ
    validate SRR_ID         Validate SRA file integrity
    check-layout SRR_ID     Check if library is single or paired-end

OPTIONS:
    -t, --threads NUM       Number of threads (default: $DEFAULT_THREADS)
    -o, --output DIR        Output directory (default: $DEFAULT_OUTPUT_DIR)
    -h, --help             Show this help message

EXAMPLES:
    $0 install
    $0 download SRX7805363
    $0 download-list SRR_Acc_List.txt
    $0 convert SRR11184871
    $0 check-layout SRR11184871

EOF
}

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Error handling
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        error_exit "Command '$1' not found. Please install it first."
    fi
}

#=============================================================================
# INSTALLATION FUNCTIONS
#=============================================================================

# Install required dependencies
install_dependencies() {
    log "Installing SRA processing dependencies..."
    
    # Check if conda is available
    if command -v conda &> /dev/null; then
        log "Installing via conda..."
        conda install -c bioconda entrez-direct sra-tools parallel-fastq-dump pigz -y
    else
        log "Conda not found. Please install dependencies manually:"
        log "  - entrez-direct"
        log "  - sra-tools"
        log "  - parallel-fastq-dump"
        log "  - pigz"
        exit 1
    fi
    
    log "Dependencies installed successfully."
}

#=============================================================================
# SRA QUERY AND DOWNLOAD FUNCTIONS
#=============================================================================

# Convert SRX ID to SRR ID(s)
srx_to_srr() {
    local srx_id="$1"
    
    log "Converting SRX ID '$srx_id' to SRR ID(s)..."
    
    check_command "esearch"
    check_command "efetch"
    
    local srr_ids
    srr_ids=$(esearch -db sra -query "$srx_id" | efetch -format runinfo | cut -d ',' -f 1 | grep SRR)
    
    if [[ -z "$srr_ids" ]]; then
        error_exit "No SRR IDs found for SRX ID: $srx_id"
    fi
    
    echo "$srr_ids"
}

# Download SRA file using prefetch
download_sra() {
    local srr_id="$1"
    
    log "Downloading SRA file for: $srr_id"
    
    check_command "prefetch"
    
    if time prefetch "$srr_id"; then
        log "Successfully downloaded: $srr_id"
    else
        error_exit "Failed to download: $srr_id"
    fi
}

# Download multiple SRA files from a list
download_sra_batch() {
    local file_list="$1"
    
    if [[ ! -f "$file_list" ]]; then
        error_exit "File list not found: $file_list"
    fi
    
    log "Downloading SRA files from list: $file_list"
    
    check_command "xargs"
    check_command "prefetch"
    
    if xargs -n1 prefetch < "$file_list"; then
        log "Batch download completed successfully."
    else
        error_exit "Batch download failed."
    fi
}

#=============================================================================
# SRA VALIDATION AND CONVERSION FUNCTIONS
#=============================================================================

# Validate SRA file integrity
validate_sra() {
    local srr_id="$1"
    
    log "Validating SRA file integrity for: $srr_id"
    
    check_command "vdb-validate"
    
    if vdb-validate "$srr_id"; then
        log "Validation successful for: $srr_id"
        return 0
    else
        log "Validation failed for: $srr_id"
        return 1
    fi
}

# Check library layout (single or paired-end)
check_library_layout() {
    local srr_id="$1"
    
    log "Checking library layout for: $srr_id"
    
    check_command "esearch"
    check_command "efetch"
    
    local layout
    layout=$(esearch -db sra -query "$srr_id" | efetch -format runinfo | cut -d ',' -f 16)
    
    if [[ -z "$layout" ]]; then
        # Alternative method
        layout=$(efetch -db sra -id "$srr_id" -format docsum | grep "LIBRARY_LAYOUT" -A 2 -m 2 | grep -v "LIBRARY_LAYOUT")
    fi
    
    log "Library layout for $srr_id: $layout"
    echo "$layout"
}

# Convert SRA to FASTQ
convert_sra_to_fastq() {
    local srr_id="$1"
    local threads="${2:-$DEFAULT_THREADS}"
    local output_dir="${3:-$DEFAULT_OUTPUT_DIR}"
    local output_prefix="${4:-$srr_id}"
    
    log "Converting SRA to FASTQ: $srr_id"
    log "Threads: $threads, Output dir: $output_dir, Prefix: $output_prefix"
    
    check_command "fasterq-dump"
    
    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"
    
    # Convert SRA to FASTQ
    if time fasterq-dump "$srr_id" \
        -e "$threads" \
        --split-files \
        -O "$output_dir" \
        -o "$output_prefix"; then
        log "Successfully converted $srr_id to FASTQ"
    else
        error_exit "Failed to convert $srr_id to FASTQ"
    fi
    
    # Compress FASTQ files
    compress_fastq_files "$output_dir"
}

# Compress FASTQ files using pigz
compress_fastq_files() {
    local directory="$1"
    local threads="${2:-8}"
    
    log "Compressing FASTQ files in: $directory"
    
    check_command "pigz"
    
    # Find and compress all .fastq files
    find "$directory" -name "*.fastq" -exec pigz -p "$threads" {} \;
    
    log "Compression completed."
}

# Batch convert multiple SRA files
batch_convert_sra() {
    local file_list="$1"
    local threads="${2:-$DEFAULT_THREADS}"
    local output_dir="${3:-$DEFAULT_OUTPUT_DIR}"
    
    if [[ ! -f "$file_list" ]]; then
        error_exit "File list not found: $file_list"
    fi
    
    log "Batch converting SRA files from list: $file_list"
    
    while IFS= read -r srr_id; do
        if [[ -n "$srr_id" && ! "$srr_id" =~ ^[[:space:]]*# ]]; then
            log "Processing: $srr_id"
            convert_sra_to_fastq "$srr_id" "$threads" "$output_dir"
        fi
    done < "$file_list"
    
    log "Batch conversion completed."
}

#=============================================================================
# PARALLEL PROCESSING FUNCTIONS
#=============================================================================

# Parallel download and conversion using GNU parallel
parallel_process_sra() {
    local file_list="$1"
    local max_jobs="${2:-3}"
    local threads="${3:-$DEFAULT_THREADS}"
    
    if [[ ! -f "$file_list" ]]; then
        error_exit "File list not found: $file_list"
    fi
    
    log "Parallel processing SRA files with $max_jobs concurrent jobs"
    
    check_command "parallel"
    
    # Using parallel-fastq-dump if available
    if command -v parallel-fastq-dump &> /dev/null; then
        while IFS= read -r srr_id; do
            if [[ -n "$srr_id" && ! "$srr_id" =~ ^[[:space:]]*# ]]; then
                parallel-fastq-dump -s "$srr_id" -t "$threads" -O "$DEFAULT_OUTPUT_DIR" --gzip
            fi
        done < "$file_list"
    else
        # Using GNU parallel with fastq-dump
        parallel -j "$max_jobs" fastq-dump --gzip --split-files {} :::: "$file_list"
    fi
}

#=============================================================================
# WORKFLOW FUNCTIONS
#=============================================================================

# Complete workflow: SRX to FASTQ
complete_workflow() {
    local srx_id="$1"
    local threads="${2:-$DEFAULT_THREADS}"
    local output_dir="${3:-$DEFAULT_OUTPUT_DIR}"
    
    log "Starting complete workflow for SRX: $srx_id"
    
    # Step 1: Convert SRX to SRR
    local srr_ids
    srr_ids=$(srx_to_srr "$srx_id")
    
    # Step 2: Process each SRR ID
    for srr_id in $srr_ids; do
        log "Processing SRR ID: $srr_id"
        
        # Download
        download_sra "$srr_id"
        
        # Validate
        if ! validate_sra "$srr_id"; then
            log "WARNING: Validation failed for $srr_id, but continuing..."
        fi
        
        # Check layout
        check_library_layout "$srr_id"
        
        # Convert to FASTQ
        convert_sra_to_fastq "$srr_id" "$threads" "$output_dir"
    done
    
    log "Complete workflow finished for SRX: $srx_id"
}

#=============================================================================
# MAIN FUNCTION
#=============================================================================

main() {
    # Parse command line arguments
    local threads="$DEFAULT_THREADS"
    local output_dir="$DEFAULT_OUTPUT_DIR"
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -o|--output)
                output_dir="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            install)
                install_dependencies
                exit 0
                ;;
            download)
                if [[ -n "${2:-}" ]]; then
                    if [[ "$2" =~ ^SRX ]]; then
                        complete_workflow "$2" "$threads" "$output_dir"
                    else
                        download_sra "$2"
                    fi
                    exit 0
                else
                    error_exit "Please provide SRX or SRR ID"
                fi
                ;;
            download-list)
                if [[ -n "${2:-}" ]]; then
                    download_sra_batch "$2"
                    exit 0
                else
                    error_exit "Please provide file list"
                fi
                ;;
            convert)
                if [[ -n "${2:-}" ]]; then
                    convert_sra_to_fastq "$2" "$threads" "$output_dir"
                    exit 0
                else
                    error_exit "Please provide SRR ID"
                fi
                ;;
            validate)
                if [[ -n "${2:-}" ]]; then
                    validate_sra "$2"
                    exit 0
                else
                    error_exit "Please provide SRR ID"
                fi
                ;;
            check-layout)
                if [[ -n "${2:-}" ]]; then
                    check_library_layout "$2"
                    exit 0
                else
                    error_exit "Please provide SRR ID"
                fi
                ;;
            *)
                error_exit "Unknown option: $1. Use -h for help."
                ;;
        esac
    done
    
    # If no command provided, show usage
    usage
    exit 1
}

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
