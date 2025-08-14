#!/usr/bin/bash

# Comment out the scripts not in use

# Install dependencies
./sra_script.sh install

# Download and process a single SRX
./sra_script.sh download SRX7805363

# Download from a list file
./sra_script.sh download-list SRR_Acc_List.txt

# Convert specific SRR to FASTQ
./sra_script.sh convert SRR11184871

# Check library layout
./sra_script.sh check-layout SRR11184871

# Use custom threads and output directory
./sra_script.sh -t 16 -o /path/to/output convert SRR11184871
