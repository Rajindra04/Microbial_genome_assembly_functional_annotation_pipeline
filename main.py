import argparse
import logging
import os
import sys
from pathlib import Path
import metagenomics_pipeline
import genomic_pipeline
import visualize
import re

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

def sanitize_organism_name(name):
    """Sanitize organism name by removing parentheses, taxon IDs, and replacing invalid characters."""
    name = re.sub(r'\(taxid_\d+\)', '', name)
    name = re.sub(r'[^a-zA-Z0-9]+', '_', name).strip('_')
    return name

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
    kraken_output_dir = os.path.join(output_dir, "kraken2_output")
    contigs = os.path.join(output_dir, "spades_output", "contigs.fasta")
    # Collect segregated FASTQ files from bin subdirectories
    fastq_pairs = []
    for bin_dir in Path(bins_dir).glob("bin*"):
        if not bin_dir.is_dir():
            continue
        for fastq1 in bin_dir.glob("*_R1.fastq"):
            organism = fastq1.stem.replace("_R1", "")
            fastq2 = fastq1.parent / f"{organism}_R2.fastq"
            if fastq2.exists():
                # Skip invalid organism names
                if organism.lower() in ['clean', 'unclassified', 'filtered', 'trimmed'] or organism.isdigit():
                    logger.warning(f"Skipping invalid organism name in {bin_dir}: {organism}")
                    continue
                fastq_pairs.append((str(fastq1), str(fastq2), organism, bin_dir.name))
                logger.info(f"Found FASTQ pair for organism {organism} in bin {bin_dir.name}: {fastq1}, {fastq2}")
    # Log bin and taxonomy information
    bin_files = list(Path(bins_dir).glob("bin*.fa"))
    logger.info(f"Found {len(bin_files)} bins in {bins_dir}: {', '.join(b.stem for b in bin_files)}")
    kraken_output_files = list(Path(kraken_output_dir).glob("*_kraken2_output.txt"))
    if kraken_output_files:
        taxa = {}
        for kraken_output in kraken_output_files:
            bin_id = Path(kraken_output).stem.replace('_kraken2_output', '')
            with open(kraken_output, 'r') as f:
                bin_taxa = set(line.strip().split('\t')[2] for line in f if line.strip().split('\t')[0] == 'C')
                taxa[bin_id] = bin_taxa
            logger.info(f"Bin {bin_id} taxa: {', '.join(bin_taxa) if bin_taxa else 'No classified taxa'}")
        logger.info("For each bin, run the genomic pipeline with an appropriate reference genome using: "
                    "python main.py --pipeline genomic --bin_id <bin_id> --ref_dir <path_to_reference> ...")
    else:
        logger.warning(f"No Kraken2 output files found in {kraken_output_dir}")
    if not fastq_pairs:
        logger.warning("No valid segregated FASTQ files found in bin subdirectories")
    return bins_dir, kraken_output_files, contigs, fastq_pairs

def run_genomic_pipeline(args, output_dir):
    """Run the genomic pipeline for a specific bin."""
    logger.info(f"Starting genomic pipeline for bin {args.bin_id}")
    bins_dir = os.path.join(output_dir, "bins")
    # Find the FASTQ pair for the specified bin
    bin_dir = Path(bins_dir) / args.bin_id
    if not bin_dir.is_dir():
        logger.error(f"Bin directory {bin_dir} does not exist")
        sys.exit(1)
    
    fastq_pairs = []
    for fastq1 in bin_dir.glob("*_R1.fastq"):
        organism = fastq1.stem.replace("_R1", "")
        fastq2 = fastq1.parent / f"{organism}_R2.fastq"
        if fastq2.exists():
            if organism.lower() in ['clean', 'unclassified', 'filtered', 'trimmed'] or organism.isdigit():
                logger.warning(f"Skipping invalid organism name in {bin_dir}: {organism}")
                continue
            fastq_pairs.append((str(fastq1), str(fastq2), organism, args.bin_id))
            logger.info(f"Found FASTQ pair for organism {organism} in bin {args.bin_id}: {fastq1}, {fastq2}")
    
    if not fastq_pairs:
        logger.error(f"No valid FASTQ pairs found for bin {args.bin_id}")
        sys.exit(1)
    
    # Ensure reference files exist
    reference_fasta = os.path.join(args.ref_dir, "reference.fasta")
    genbank_file = os.path.join(args.ref_dir, "reference.gbk")
    for f, desc in [(reference_fasta, "Reference FASTA"), (genbank_file, "GenBank")]:
        if not os.path.isfile(f):
            logger.error(f"{desc} ({f}) not found in {args.ref_dir}")
            sys.exit(1)
    
    genomic_output_dirs = []
    for fastq1, fastq2, organism, bin_id in fastq_pairs:
        logger.info(f"Processing organism {organism} from bin {bin_id}")
        sanitized_organism = sanitize_organism_name(organism)
        organism_output_dir = os.path.join(output_dir, "genomic", f"{bin_id}_{sanitized_organism}")
        os.makedirs(organism_output_dir, exist_ok=True)
        
        # Construct command-line arguments for genomic_pipeline
        sys.argv = [
            sys.argv[0],
            "-s", sanitized_organism,
            "-1", fastq1,
            "-2", fastq2,
            "-r", args.ref_dir,
            "-d", args.snpeff_db,
            "-o", organism_output_dir,
            "-t", str(args.threads),
            "-k", args.organism
        ]
        if args.keep_temp:
            sys.argv.append("--keep-temp")
        if args.dry_run:
            sys.argv.append("--dry-run")
        if args.verbose:
            sys.argv.append("--verbose")
        if args.resume:
            sys.argv.append("--resume")
        if args.skip_emapper:
            sys.argv.append("--skip-emapper")
        
        try:
            genomic_pipeline.main()
            genomic_output_dirs.append(organism_output_dir)
        except Exception as e:
            logger.error(f"Genomic pipeline failed for organism {organism} in bin {bin_id}: {str(e)}")
            continue
    
    logger.info(f"Genomic pipeline completed for {len(genomic_output_dirs)} organisms in bin {args.bin_id}")
    return genomic_output_dirs

def run_visualization(args, output_dir):
    """Run the visualization pipeline for genomic outputs."""
    logger.info("Starting visualization pipeline")
    genomic_output_dirs = []
    
    if args.bin_id:
        # Process only the specified bin
        genomic_output_dirs = [d for d in Path(output_dir).glob(f"genomic/{args.bin_id}_*") if d.is_dir()]
        if not genomic_output_dirs:
            logger.error(f"No genomic output directory found for bin {args.bin_id}")
            sys.exit(1)
    else:
        # Process all genomic outputs
        genomic_output_dirs = [d for d in Path(output_dir).glob("genomic/*") if d.is_dir()]
        if not genomic_output_dirs:
            logger.warning("No genomic output directories found for visualization")
            return
    
    for organism_output_dir in genomic_output_dirs:
        # Extract the sanitized organism name from the directory name
        organism = Path(organism_output_dir).name.split('_', 1)[1] if '_' in Path(organism_output_dir).name else Path(organism_output_dir).name
        annot_file = os.path.join(organism_output_dir, f"{organism}.annotations.txt")
        if not os.path.exists(annot_file):
            logger.warning(f"Annotation file {annot_file} not found, skipping visualization")
            continue
        
        logger.info(f"Generating visualizations for {organism}")
        vis_output_dir = os.path.join(output_dir, "visualizations", organism)
        os.makedirs(vis_output_dir, exist_ok=True)
        
        # Update sys.argv for visualize.py to include --output-dir
        sys.argv = [sys.argv[0], annot_file, "--output-dir", vis_output_dir]
        try:
            visualize.main()
        except Exception as e:
            logger.error(f"Visualization failed for {organism}: {str(e)}")
            continue
    
    logger.info("Visualization pipeline completed")

def main():
    parser = argparse.ArgumentParser(description="Integrated metagenomics and genomic analysis pipeline")
    parser.add_argument("--pipeline", choices=["metagenomics", "genomic", "visualize"], required=True,
                        help="Pipeline to run: metagenomics, genomic, or visualize")
    parser.add_argument("--fastq1", help="Path to forward FASTQ file (required for metagenomics)")
    parser.add_argument("--fastq2", help="Path to reverse FASTQ file (required for metagenomics)")
    parser.add_argument("--output", default="output", help="Output directory")
    parser.add_argument("--ref_dir", help="Directory containing reference.fasta and reference.gbk (required for genomic)")
    parser.add_argument("--snpeff_db", help="SnpEff database name (required for genomic)")
    parser.add_argument("--db_path", help="Path to Kraken2 database (required for metagenomics if KRAKEN2_DB_PATH not set)")
    parser.add_argument("--bin_id", help="Bin ID for genomic or visualize pipeline (e.g., bin1)")
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

    if args.pipeline == "metagenomics":
        if not (args.fastq1 and args.fastq2):
            logger.error("Both --fastq1 and --fastq2 are required for metagenomics pipeline")
            sys.exit(1)
        if not args.db_path and not os.environ.get("KRAKEN2_DB_PATH"):
            logger.error(
                "Kraken2 database path must be specified via --db_path or the KRAKEN2_DB_PATH environment variable "
                "for metagenomics pipeline. Set it with 'export KRAKEN2_DB_PATH=/path/to/database' or provide --db_path."
            )
            sys.exit(1)
        if args.db_path and not os.path.exists(args.db_path):
            logger.error(f"Kraken2 database path {args.db_path} does not exist.")
            sys.exit(1)
        run_metagenomics_pipeline(args, output_dir)
    elif args.pipeline == "genomic":
        if not (args.ref_dir and args.snpeff_db):
            logger.error("Both --ref_dir and --snpeff_db are required for genomic pipeline")
            sys.exit(1)
        if not args.bin_id:
            logger.error("--bin_id is required for genomic pipeline to specify which bin to process")
            sys.exit(1)
        run_genomic_pipeline(args, output_dir)
    elif args.pipeline == "visualize":
        run_visualization(args, output_dir)

if __name__ == "__main__":
    main()

