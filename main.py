import argparse
import logging
import os
import sys
from pathlib import Path
import metagenomics_pipeline
import genomic_pipeline
import visualize

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pipeline.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_metagenomics_pipeline(args, output_dir):
    """Run the metagenomics pipeline."""
    logger.info("Starting metagenomics pipeline")
    metagenomics_pipeline.main(
        input_fastq1=args.fastq1,
        input_fastq2=args.fastq2,
        output_dir=output_dir,
        resume=args.resume,
        db_path=args.db_path
    )
    logger.info("Metagenomics pipeline completed")
    # Return paths to segregated FASTQ files and contigs
    bins_dir = os.path.join(output_dir, "bins")
    kraken_output = os.path.join(output_dir, "kraken2_output", "kraken2_output.txt")
    contigs = os.path.join(output_dir, "spades_output", "contigs.fasta")
    # Verify segregated FASTQ files exist
    fastq_pairs = []
    for fastq1 in Path(output_dir).glob("*_R1.fastq"):
        organism = fastq1.stem.replace("_R1", "")
        fastq2 = fastq1.parent / f"{organism}_R2.fastq"
        if fastq2.exists():
            fastq_pairs.append((str(fastq1), str(fastq2), organism))
    if not fastq_pairs:
        logger.warning("No segregated FASTQ files found in output directory")
    return bins_dir, kraken_output, contigs, fastq_pairs

def run_genomic_pipeline(args, fastq_pairs, output_dir):
    """Run the genomic pipeline for each bin's FASTQ files."""
    logger.info("Starting genomic pipeline")
    genomic_output_dirs = []
    
    if not fastq_pairs:
        logger.error("No segregated FASTQ pairs available for genomic pipeline")
        return []
    
    # Ensure reference files exist
    reference_fasta = os.path.join(args.ref_dir, "reference.fasta")
    genbank_file = os.path.join(args.ref_dir, "reference.gbk")
    for f, desc in [(reference_fasta, "Reference FASTA"), (genbank_file, "GenBank")]:
        if not os.path.isfile(f):
            logger.error(f"{desc} ({f}) not found.")
            sys.exit(1)
    
    for fastq1, fastq2, organism in fastq_pairs:
        logger.info(f"Processing organism: {organism}")
        organism_output_dir = os.path.join(output_dir, "genomic", organism)
        os.makedirs(organism_output_dir, exist_ok=True)
        
        # Prepare arguments for genomic pipeline
        genomic_args = argparse.Namespace(
            sample=organism,
            fastq1=fastq1,
            fastq2=fastq2,
            ref_dir=args.ref_dir,
            snpeff_db=args.snpeff_db,
            output_dir=organism_output_dir,
            threads=args.threads,
            organism=args.organism,
            keep_temp=args.keep_temp,
            dry_run=args.dry_run,
            verbose=args.verbose,
            resume=args.resume,
            skip_emapper=args.skip_emapper
        )
        
        try:
            genomic_pipeline.main(genomic_args)
            genomic_output_dirs.append(organism_output_dir)
        except Exception as e:
            logger.error(f"Genomic pipeline failed for organism {organism}: {str(e)}")
            continue
    
    logger.info(f"Genomic pipeline completed for {len(genomic_output_dirs)} organisms")
    return genomic_output_dirs

def run_visualization(args, genomic_output_dirs, output_dir):
    """Run the visualization pipeline for each genomic output."""
    logger.info("Starting visualization pipeline")
    
    for organism_output_dir in genomic_output_dirs:
        annot_file = os.path.join(organism_output_dir, f"{Path(organism_output_dir).name}.annotations.txt")
        if not os.path.exists(annot_file):
            logger.warning(f"Annotation file {annot_file} not found, skipping visualization")
            continue
        
        logger.info(f"Generating visualizations for {Path(organism_output_dir).name}")
        vis_output_dir = os.path.join(output_dir, "visualizations", Path(organism_output_dir).name)
        os.makedirs(vis_output_dir, exist_ok=True)
        
        # Update sys.argv for visualize.py
        sys.argv = [sys.argv[0], annot_file]
        try:
            visualize.main()
        except Exception as e:
            logger.error(f"Visualization failed for {Path(organism_output_dir).name}: {str(e)}")
            continue
    
    logger.info("Visualization pipeline completed")

def main():
    parser = argparse.ArgumentParser(description="Integrated metagenomics and genomic analysis pipeline")
    parser.add_argument("--pipeline", choices=["all", "metagenomics", "genomic", "visualize"], default="all",
                        help="Pipeline to run: all, metagenomics, genomic, or visualize")
    parser.add_argument("--fastq1", help="Path to forward FASTQ file (required for metagenomics/all)")
    parser.add_argument("--fastq2", help="Path to reverse FASTQ file (required for metagenomics/all)")
    parser.add_argument("--output", default="output", help="Output directory")
    parser.add_argument("--ref_dir", help="Directory containing reference.fasta and reference.gbk (required for genomic/all)")
    parser.add_argument("--snpeff_db", help="SnpEff database name (required for genomic/all)")
    parser.add_argument("--db_path", help="Path to Kraken2 database (overrides KRAKEN2_DB_PATH)")
    parser.add_argument("--threads", type=int, default=8, help="Number of CPU threads")
    parser.add_argument("--organism", default="eco", help="Organism code for KEGG")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--resume", action="store_true", help="Resume from completed steps")
    parser.add_argument("--skip-emapper", action="store_true", help="Skip eggNOG-mapper step if database is unavailable")
    args = parser.parse_args()

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    fastq_pairs = []
    if args.pipeline in ["all", "metagenomics"]:
        if not (args.fastq1 and args.fastq2):
            logger.error("Both --fastq1 and --fastq2 are required for metagenomics pipeline")
            sys.exit(1)
        bins_dir, kraken_output, contigs, fastq_pairs = run_metagenomics_pipeline(args, output_dir)
    else:
        bins_dir, contigs = None, None

    if args.pipeline in ["all", "genomic"]:
        if not (args.ref_dir and args.snpeff_db):
            logger.error("Both --ref_dir and --snpeff_db are required for genomic pipeline")
            sys.exit(1)
        if not fastq_pairs and args.pipeline == "genomic":
            # Look for existing segregated FASTQ files if running genomic pipeline standalone
            fastq_pairs = []
            for fastq1 in Path(output_dir).glob("*_R1.fastq"):
                organism = fastq1.stem.replace("_R1", "")
                fastq2 = fastq1.parent / f"{organism}_R2.fastq"
                if fastq2.exists():
                    fastq_pairs.append((str(fastq1), str(fastq2), organism))
            if not fastq_pairs:
                logger.error("No segregated FASTQ files found for genomic pipeline")
                sys.exit(1)
        genomic_output_dirs = run_genomic_pipeline(args, fastq_pairs, output_dir)
    else:
        genomic_output_dirs = []

    if args.pipeline in ["all", "visualize"]:
        if not genomic_output_dirs:
            # If genomic pipeline wasn't run, look for existing genomic outputs
            genomic_output_dirs = [d for d in Path(output_dir).glob("genomic/*") if d.is_dir()]
        if not genomic_output_dirs:
            logger.warning("No genomic output directories found for visualization")
        else:
            run_visualization(args, genomic_output_dirs, output_dir)

if __name__ == "__main__":
    main()