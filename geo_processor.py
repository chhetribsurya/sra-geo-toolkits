#!/usr/bin/env python3

"""
GEO Dataset Download and Processing Script (Python)
Dependencies: GEOparse, pandas, numpy
"""

import os
import sys
import argparse
import logging
from typing import List, Dict, Optional, Union
from pathlib import Path

import pandas as pd
import numpy as np
import GEOparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('geo_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

#=============================================================================
# GEO DATA DOWNLOAD AND PROCESSING
#=============================================================================

class GEODataProcessor:
    """
    A class for downloading and processing GEO datasets.
    
    Attributes:
        output_dir (str): Base output directory for all analyses
    """
    
    def __init__(self, output_dir: str = "./geo_analysis"):
        """
        Initialize the GEO data processor.
        
        Args:
            output_dir (str): Base output directory
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"GEO Data Processor initialized with output dir: {self.output_dir}")

    def download_geo_dataset(self, gse_id: str, destdir: Optional[str] = None) -> GEOparse.GEOTypes.GSE:
        """
        Download GEO dataset using GEOparse.
        
        Args:
            gse_id (str): GSE accession number (e.g., 'GSE188486')
            destdir (str, optional): Destination directory for download
            
        Returns:
            GEOparse.GEOTypes.GSE: Downloaded GSE object
            
        Raises:
            Exception: If download fails
        """
        if destdir is None:
            destdir = str(self.output_dir / gse_id)
        
        logger.info(f"Downloading GEO dataset: {gse_id}")
        
        try:
            # Download dataset
            gse = GEOparse.get_GEO(geo=gse_id, destdir=destdir)
            logger.info(f"Successfully downloaded {gse_id}")
            
            # Print basic information
            logger.info(f"Dataset title: {gse.metadata.get('title', ['N/A'])[0]}")
            logger.info(f"Number of samples: {len(gse.phenotype_data)}")
            logger.info(f"Number of supplementary files: {len(gse.supplementary_files)}")
            
            return gse
            
        except Exception as e:
            logger.error(f"Error downloading {gse_id}: {str(e)}")
            raise

    def download_supplementary_files(self, gse: GEOparse.GEOTypes.GSE, 
                                   output_dir: Optional[str] = None) -> List[str]:
        """
        Download all supplementary files for a GSE dataset.
        
        Args:
            gse (GEOparse.GEOTypes.GSE): GSE object
            output_dir (str, optional): Output directory for files
            
        Returns:
            List[str]: List of downloaded file names
        """
        if output_dir is None:
            output_dir = str(self.output_dir)
        
        downloaded_files = []
        
        if not gse.supplementary_files:
            logger.info("No supplementary files available for download")
            return downloaded_files
        
        logger.info(f"Downloading {len(gse.supplementary_files)} supplementary files...")
        
        for supp_file in gse.supplementary_files:
            try:
                filename = supp_file.split('/')[-1]
                local_path = os.path.join(output_dir, filename)
                
                logger.info(f"Downloading: {filename}")
                GEOparse.utils.download_from_url(supp_file, filename=local_path)
                downloaded_files.append(filename)
                logger.info(f"Successfully downloaded: {filename}")
                
            except Exception as e:
                logger.error(f"Failed to download {supp_file}: {str(e)}")
        
        logger.info(f"Downloaded {len(downloaded_files)} supplementary files")
        return downloaded_files

    def extract_metadata(self, gse: GEOparse.GEOTypes.GSE, 
                        selected_columns: Optional[List[str]] = None,
                        output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Extract and process metadata from GSE object.
        
        Args:
            gse (GEOparse.GEOTypes.GSE): GSE object
            selected_columns (List[str], optional): Specific columns to extract
            output_file (str, optional): Output file path
            
        Returns:
            pd.DataFrame: Processed metadata
        """
        logger.info("Extracting metadata...")
        
        # Get metadata as DataFrame
        metadata_df = gse.phenotype_data.copy()
        
        logger.info(f"Total metadata columns: {len(metadata_df.columns)}")
        logger.info(f"Available columns: {list(metadata_df.columns)}")
        
        # Select specific columns if requested
        if selected_columns:
            available_cols = [col for col in selected_columns if col in metadata_df.columns]
            missing_cols = [col for col in selected_columns if col not in metadata_df.columns]
            
            if available_cols:
                metadata_df = metadata_df[available_cols]
                logger.info(f"Selected {len(available_cols)} columns from metadata: {available_cols}")
            else:
                logger.warning("None of the specified columns found in metadata")
            
            if missing_cols:
                logger.warning(f"Missing columns: {missing_cols}")
        
        # Save to file if specified
        if output_file:
            metadata_df.to_csv(output_file, sep="\t", index=True)
            logger.info(f"Metadata saved to: {output_file}")
        
        return metadata_df

    def filter_samples_by_criteria(self, metadata_df: pd.DataFrame, 
                                 filter_column: str, 
                                 filter_pattern: str,
                                 case_sensitive: bool = False) -> pd.DataFrame:
        """
        Filter samples based on metadata criteria.
        
        Args:
            metadata_df (pd.DataFrame): Metadata DataFrame
            filter_column (str): Column name to filter on
            filter_pattern (str): Pattern to match (supports regex)
            case_sensitive (bool): Whether filtering is case sensitive
            
        Returns:
            pd.DataFrame: Filtered metadata
            
        Raises:
            ValueError: If filter column not found
        """
        if filter_column not in metadata_df.columns:
            available_cols = list(metadata_df.columns)
            raise ValueError(f"Column '{filter_column}' not found in metadata. Available columns: {available_cols}")
        
        logger.info(f"Filtering samples on '{filter_column}' with pattern: '{filter_pattern}'")
        
        # Apply filter
        if case_sensitive:
            mask = metadata_df[filter_column].str.contains(filter_pattern, na=False, regex=True)
        else:
            mask = metadata_df[filter_column].str.contains(filter_pattern, case=False, na=False, regex=True)
        
        filtered_df = metadata_df[mask].copy()
        
        logger.info(f"Filtered from {len(metadata_df)} to {len(filtered_df)} samples")
        
        if len(filtered_df) == 0:
            logger.warning(f"No samples found matching pattern '{filter_pattern}' in column '{filter_column}'")
        
        return filtered_df

    def get_unique_values(self, metadata_df: pd.DataFrame, column: str) -> pd.Series:
        """
        Get unique values from a metadata column.
        
        Args:
            metadata_df (pd.DataFrame): Metadata DataFrame
            column (str): Column name to analyze
            
        Returns:
            pd.Series: Unique values and their counts
            
        Raises:
            ValueError: If column not found
        """
        if column not in metadata_df.columns:
            raise ValueError(f"Column '{column}' not found in metadata")
        
        unique_values = metadata_df[column].value_counts()
        logger.info(f"Unique values in '{column}' column:")
        for value, count in unique_values.items():
            logger.info(f"  {value}: {count}")
        
        return unique_values

    def create_sample_annotation(self, metadata_df: pd.DataFrame,
                               key_columns: Optional[List[str]] = None,
                               output_file: Optional[str] = None) -> pd.DataFrame:
        """
        Create a cleaned sample annotation file from metadata.
        
        Args:
            metadata_df (pd.DataFrame): Metadata DataFrame
            key_columns (List[str], optional): Important columns to include
            output_file (str, optional): Output file path
            
        Returns:
            pd.DataFrame: Sample annotation DataFrame
        """
        logger.info("Creating sample annotation...")
        
        # Default key columns if not specified
        if key_columns is None:
            key_columns = [
                "title", "geo_accession", "source_name_ch1", 
                "characteristics_ch1", "description", "supplementary_file"
            ]
        
        # Select available key columns
        available_cols = [col for col in key_columns if col in metadata_df.columns]
        
        if available_cols:
            annotation_df = metadata_df[available_cols].copy()
            logger.info(f"Created annotation with {len(available_cols)} columns: {available_cols}")
        else:
            annotation_df = metadata_df.copy()
            logger.info("Using all metadata columns for annotation")
        
        # Clean up column names (replace spaces with underscores)
        annotation_df.columns = [col.replace(' ', '_') for col in annotation_df.columns]
        
        # Save to file if specified
        if output_file:
            annotation_df.to_csv(output_file, sep="\t", index=True)
            logger.info(f"Sample annotation saved to: {output_file}")
        
        return annotation_df

    def rename_supplementary_files(self, gse: GEOparse.GEOTypes.GSE,
                                 metadata_df: pd.DataFrame,
                                 source_dir: str,
                                 output_dir: str,
                                 identifier_column: str = "source_name_ch1") -> Dict[str, str]:
        """
        Rename supplementary files based on sample metadata.
        
        Args:
            gse (GEOparse.GEOTypes.GSE): GSE object
            metadata_df (pd.DataFrame): Metadata DataFrame
            source_dir (str): Directory containing downloaded files
            output_dir (str): Directory for renamed files
            identifier_column (str): Column to use for renaming
            
        Returns:
            Dict[str, str]: Mapping of old filenames to new filenames
        """
        logger.info("Renaming supplementary files based on metadata...")
        
        source_path = Path(source_dir)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        renamed_files = {}
        
        # Create mapping from supplementary files to sample information
        for index, row in metadata_df.iterrows():
            try:
                geo_accession = row['geo_accession']
                
                # Get identifier for renaming
                if identifier_column in row and pd.notna(row[identifier_column]):
                    identifier = str(row[identifier_column]).replace(' ', '_').replace('/', '_')
                else:
                    identifier = geo_accession
                
                # Find supplementary file URL for this sample
                supp_file_url = None
                for col in row.index:
                    if 'supplementary_file' in col and pd.notna(row[col]):
                        supp_file_url = row[col]
                        break
                
                if supp_file_url:
                    # Extract original filename
                    original_filename = supp_file_url.split('/')[-1]
                    original_path = source_path / original_filename
                    
                    if original_path.exists():
                        # Create new filename
                        file_extension = original_path.suffix
                        new_filename = f"{identifier}-{original_filename}"
                        new_path = output_path / new_filename
                        
                        # Copy and rename file
                        import shutil
                        shutil.copy2(original_path, new_path)
                        renamed_files[original_filename] = new_filename
                        
                        logger.info(f"Renamed: {original_filename} -> {new_filename}")
                    else:
                        logger.warning(f"Original file not found: {original_filename}")
                        
            except Exception as e:
                logger.error(f"Error renaming file for {index}: {str(e)}")
        
        logger.info(f"Renamed {len(renamed_files)} files")
        return renamed_files

    def generate_dataset_summary(self, gse: GEOparse.GEOTypes.GSE,
                               metadata_df: pd.DataFrame,
                               output_file: Optional[str] = None) -> Dict:
        """
        Generate a comprehensive summary of the dataset.
        
        Args:
            gse (GEOparse.GEOTypes.GSE): GSE object
            metadata_df (pd.DataFrame): Metadata DataFrame
            output_file (str, optional): Output file for summary
            
        Returns:
            Dict: Dataset summary information
        """
        logger.info("Generating dataset summary...")
        
        # Basic dataset information
        summary = {
            'gse_id': gse.name,
            'title': gse.metadata.get('title', ['N/A'])[0],
            'summary': gse.metadata.get('summary', ['N/A'])[0],
            'overall_design': gse.metadata.get('overall_design', ['N/A'])[0],
            'submission_date': gse.metadata.get('submission_date', ['N/A'])[0],
            'last_update_date': gse.metadata.get('last_update_date', ['N/A'])[0],
            'pubmed_id': gse.metadata.get('pubmed_id', ['N/A'])[0],
            'platform_count': len(gse.gpls),
            'sample_count': len(metadata_df),
            'supplementary_file_count': len(gse.supplementary_files),
            'columns_in_metadata': list(metadata_df.columns)
        }
        
        # Platform information
        if gse.gpls:
            summary['platforms'] = {}
            for platform_id, platform in gse.gpls.items():
                summary['platforms'][platform_id] = {
                    'title': platform.metadata.get('title', ['N/A'])[0],
                    'technology': platform.metadata.get('technology', ['N/A'])[0],
                    'organism': platform.metadata.get('organism', ['N/A'])[0]
                }
        
        # Sample type distribution
        if 'source_name_ch1' in metadata_df.columns:
            sample_types = metadata_df['source_name_ch1'].value_counts().to_dict()
            summary['sample_type_distribution'] = sample_types
        
        # Print summary
        logger.info("=== Dataset Summary ===")
        logger.info(f"GSE ID: {summary['gse_id']}")
        logger.info(f"Title: {summary['title']}")
        logger.info(f"Samples: {summary['sample_count']}")
        logger.info(f"Platforms: {summary['platform_count']}")
        logger.info(f"Supplementary files: {summary['supplementary_file_count']}")
        
        if 'sample_type_distribution' in summary:
            logger.info("Sample type distribution:")
            for sample_type, count in summary['sample_type_distribution'].items():
                logger.info(f"  {sample_type}: {count}")
        
        # Save summary to file if specified
        if output_file:
            import json
            with open(output_file, 'w') as f:
                json.dump(summary, f, indent=2, default=str)
            logger.info(f"Summary saved to: {output_file}")
        
        return summary

#=============================================================================
# BATCH PROCESSING
#=============================================================================

def process_multiple_datasets(gse_ids: List[str], 
                            output_base_dir: str = "./batch_geo_analysis",
                            download_supplementary: bool = True,
                            filter_pattern: Optional[str] = None,
                            filter_column: str = "title") -> Dict[str, Dict]:
    """
    Process multiple GEO datasets in batch.
    
    Args:
        gse_ids (List[str]): List of GSE accession numbers
        output_base_dir (str): Base output directory
        download_supplementary (bool): Whether to download supplementary files
        filter_pattern (str, optional): Pattern to filter samples
        filter_column (str): Column to use for filtering
        
    Returns:
        Dict[str, Dict]: Results for each dataset
    """
    logger.info(f"Processing {len(gse_ids)} datasets in batch...")
    
    base_path = Path(output_base_dir)
    base_path.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    for gse_id in gse_ids:
        logger.info(f"\n{'='*50}")
        logger.info(f"Processing dataset: {gse_id}")
        logger.info(f"{'='*50}")
        
        try:
            # Create individual output directory
            dataset_dir = base_path / gse_id
            
            # Initialize processor
            processor = GEODataProcessor(str(dataset_dir))
            
            # Download dataset
            gse = processor.download_geo_dataset(gse_id)
            
            # Extract metadata
            metadata_file = dataset_dir / f"{gse_id}_metadata.tsv"
            metadata = processor.extract_metadata(gse, output_file=str(metadata_file))
            
            # Download supplementary files if requested
            supplementary_files = []
            if download_supplementary:
                supplementary_files = processor.download_supplementary_files(gse, str(dataset_dir))
            
            # Filter samples if pattern provided
            filtered_metadata = None
            if filter_pattern and filter_column in metadata.columns:
                try:
                    filtered_metadata = processor.filter_samples_by_criteria(
                        metadata, filter_column, filter_pattern
                    )
                    if len(filtered_metadata) > 0:
                        filtered_file = dataset_dir / f"{gse_id}_filtered_metadata.tsv"
                        filtered_metadata.to_csv(filtered_file, sep="\t", index=True)
                except Exception as e:
                    logger.error(f"Error filtering samples: {str(e)}")
            
            # Generate summary
            summary_file = dataset_dir / f"{gse_id}_summary.json"
            summary = processor.generate_dataset_summary(gse, metadata, str(summary_file))
            
            # Create sample annotation
            annotation_file = dataset_dir / f"{gse_id}_sample_annotation.tsv"
            annotation = processor.create_sample_annotation(metadata, output_file=str(annotation_file))
            
            # Store results
            results[gse_id] = {
                'status': 'success',
                'gse_object': gse,
                'metadata': metadata,
                'filtered_metadata': filtered_metadata,
                'supplementary_files': supplementary_files,
                'summary': summary,
                'output_directory': str(dataset_dir)
            }
            
            logger.info(f"Successfully processed {gse_id}")
            
        except Exception as e:
            logger.error(f"Error processing {gse_id}: {str(e)}")
            results[gse_id] = {
                'status': 'error',
                'error': str(e)
            }
    
    # Generate batch summary
    successful = sum(1 for r in results.values() if r['status'] == 'success')
    failed = len(results) - successful
    
    logger.info(f"\n{'='*50}")
    logger.info(f"Batch processing completed!")
    logger.info(f"Successful: {successful}/{len(gse_ids)}")
    logger.info(f"Failed: {failed}/{len(gse_ids)}")
    logger.info(f"Results saved in: {output_base_dir}")
    logger.info(f"{'='*50}")
    
    return results

#=============================================================================
# MAIN WORKFLOW FUNCTION
#=============================================================================

def analyze_geo_dataset(gse_id: str,
                       output_dir: str = "./geo_analysis",
                       download_supplementary: bool = True,
                       filter_pattern: Optional[str] = None,
                       filter_column: str = "title",
                       selected_columns: Optional[List[str]] = None) -> Dict:
    """
    Complete workflow for analyzing a single GEO dataset.
    
    Args:
        gse_id (str): GSE accession number
        output_dir (str): Output directory
        download_supplementary (bool): Whether to download supplementary files
        filter_pattern (str, optional): Pattern to filter samples
        filter_column (str): Column to use for filtering
        selected_columns (List[str], optional): Specific metadata columns to extract
        
    Returns:
        Dict: Analysis results
    """
    logger.info(f"Starting complete analysis for {gse_id}")
    
    # Initialize processor
    processor = GEODataProcessor(output_dir)
    
    # Download dataset
    gse = processor.download_geo_dataset(gse_id)
    
    # Extract metadata
    metadata_file = Path(output_dir) / f"{gse_id}_metadata.tsv"
    metadata = processor.extract_metadata(
        gse, 
        selected_columns=selected_columns,
        output_file=str(metadata_file)
    )
    
    # Download supplementary files
    supplementary_files = []
    if download_supplementary:
        supplementary_files = processor.download_supplementary_files(gse, output_dir)
    
    # Filter samples if requested
    filtered_metadata = None
    if filter_pattern:
        try:
            filtered_metadata = processor.filter_samples_by_criteria(
                metadata, filter_column, filter_pattern
            )
            if len(filtered_metadata) > 0:
                filtered_file = Path(output_dir) / f"{gse_id}_filtered_metadata.tsv"
                filtered_metadata.to_csv(filtered_file, sep="\t", index=True)
        except Exception as e:
            logger.error(f"Error filtering samples: {str(e)}")
    
    # Generate summary
    summary_file = Path(output_dir) / f"{gse_id}_summary.json"
    summary = processor.generate_dataset_summary(gse, metadata, str(summary_file))
    
    # Create sample annotation
    annotation_file = Path(output_dir) / f"{gse_id}_sample_annotation.tsv"
    annotation = processor.create_sample_annotation(metadata, output_file=str(annotation_file))
    
    results = {
        'gse_id': gse_id,
        'gse_object': gse,
        'metadata': metadata,
        'filtered_metadata': filtered_metadata,
        'supplementary_files': supplementary_files,
        'summary': summary,
        'annotation': annotation,
        'output_directory': output_dir
    }
    
    logger.info(f"Analysis completed for {gse_id}")
    return results

#=============================================================================
# COMMAND LINE INTERFACE
#=============================================================================

def main():
    """
    Main function for command line interface.
    """
    parser = argparse.ArgumentParser(
        description='GEO Dataset Download and Processing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download and analyze single dataset
  python geo_processor.py analyze GSE188486 --output-dir ./results

  # Download with filtering
  python geo_processor.py analyze GSE188486 --filter-pattern "H3K" --filter-column "title"

  # Batch processing
  python geo_processor.py batch GSE188486 GSE123456 --output-dir ./batch_results

  # Download only (no processing)
  python geo_processor.py download GSE188486 --supplementary
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Analyze command (complete workflow)
    analyze_parser = subparsers.add_parser('analyze', help='Complete analysis workflow')
    analyze_parser.add_argument('gse_id', help='GSE accession number')
    analyze_parser.add_argument('--output-dir', default='./geo_analysis', 
                               help='Output directory')
    analyze_parser.add_argument('--no-supplementary', action='store_true',
                               help='Skip downloading supplementary files')
    analyze_parser.add_argument('--filter-pattern', 
                               help='Pattern to filter samples (regex supported)')
    analyze_parser.add_argument('--filter-column', default='title',
                               help='Column to use for filtering (default: title)')
    analyze_parser.add_argument('--selected-columns', nargs='+',
                               help='Specific metadata columns to extract')
    
    # Download command (download only)
    download_parser = subparsers.add_parser('download', help='Download dataset only')
    download_parser.add_argument('gse_id', help='GSE accession number')
    download_parser.add_argument('--output-dir', default='./geo_download',
                               help='Output directory')
    download_parser.add_argument('--supplementary', action='store_true',
                               help='Download supplementary files')
    
    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Process multiple datasets')
    batch_parser.add_argument('gse_ids', nargs='+', help='GSE accession numbers')
    batch_parser.add_argument('--output-dir', default='./batch_geo_analysis',
                             help='Output directory')
    batch_parser.add_argument('--no-supplementary', action='store_true',
                             help='Skip downloading supplementary files')
    batch_parser.add_argument('--filter-pattern',
                             help='Pattern to filter samples')
    batch_parser.add_argument('--filter-column', default='title',
                             help='Column to use for filtering')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        results = analyze_geo_dataset(
            gse_id=args.gse_id,
            output_dir=args.output_dir,
            download_supplementary=not args.no_supplementary,
            filter_pattern=args.filter_pattern,
            filter_column=args.filter_column,
            selected_columns=args.selected_columns
        )
        logger.info("Analysis completed successfully!")
    
    elif args.command == 'download':
        processor = GEODataProcessor(args.output_dir)
        gse = processor.download_geo_dataset(args.gse_id)
        
        if args.supplementary:
            processor.download_supplementary_files(gse, args.output_dir)
        
        logger.info("Download completed successfully!")
    
    elif args.command == 'batch':
        results = process_multiple_datasets(
            gse_ids=args.gse_ids,
            output_base_dir=args.output_dir,
            download_supplementary=not args.no_supplementary,
            filter_pattern=args.filter_pattern,
            filter_column=args.filter_column
        )
        logger.info("Batch processing completed successfully!")
    
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
