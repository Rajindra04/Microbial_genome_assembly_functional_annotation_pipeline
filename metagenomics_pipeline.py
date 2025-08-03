import subprocess
import os
import logging
from pathlib import Path
from collections import Counter
from Bio import SeqIO
import re

# Set up logging to both file and console
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# File handler for logging to file
file_handler = logging.FileHandler('metagenomics_pipeline.log')
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(file_handler)

# Console handler for logging to screen
console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(console_handler)

def sanitize_organism_name(name):
    """Sanitize organism name by removing parentheses, taxon IDs, and replacing invalid characters."""
    # Remove taxon ID in parentheses (e.g., "(taxid_1590)")
    name = re.sub(r'\(taxid_\d+\)', '', name)
    # Replace spaces and special characters with underscores, remove trailing underscores
    name = re.sub(r'[^a-zA-Z0-9]+', '_', name).strip('_')
    return name

def run_command(command, step_name):
    """Run a shell command and handle errors."""
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, capture_output=True)
        logging.info(f"{step_name} completed successfully")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in {step_name}: {e.stderr}")
        raise

def quality_control(input_fastq1, input_fastq2, output_dir, resume=False):
    """Perform quality control using fastp."""
    output_fastq1 = os.path.join(output_dir, "clean_R1.fastq")
    output_fastq2 = os.path.join(output_dir, "clean_R2.fastq")
    
    if resume and os.path.exists(output_fastq1) and os.path.exists(output_fastq2):
        logging.info("Quality control outputs found, skipping fastp")
        return output_fastq1, output_fastq2
    
    html_report = os.path.join(output_dir, "fastp.html")
    json_report = os.path.join(output_dir, "fastp.json")
    
    command = (
        f"fastp --in1 {input_fastq1} --in2 {input_fastq2} "
        f"--out1 {output_fastq1} --out2 {output_fastq2} "
        f"--html {html_report} --json {json_report} "
        f"--thread 8"
    )
    run_command(command, "Quality control with fastp")
    return output_fastq1, output_fastq2

def denovo_assembly(clean_fastq1, clean_fastq2, output_dir, resume=False):
    """Perform de novo assembly using SPAdes."""
    spades_output = os.path.join(output_dir, "spades_output")
    contigs = os.path.join(spades_output, "contigs.fasta")
    
    if resume and os.path.exists(contigs):
        logging.info("SPAdes contigs found, skipping de novo assembly")
        return contigs
    
    command = (
        f"spades.py --meta -1 {clean_fastq1} -2 {clean_fastq2} "
        f"-o {spades_output} -t 16"
    )
    run_command(command, "De novo assembly with SPAdes")
    logging.info(f"Contigs file generated: {contigs}")
    return contigs

def bin_contigs(contigs, clean_fastq1, clean_fastq2, output_dir, resume=False):
    """Bin contigs using MetaBAT2."""
    logging.info(f"Using contigs file for binning: {contigs}")
    bam_file = os.path.join(output_dir, "aligned.bam")
    bins_dir = os.path.join(output_dir, "bins")
    
    if resume and os.path.exists(bins_dir) and any(Path(bins_dir).glob("bin*.fa")):
        logging.info("MetaBAT2 bins found, skipping binning")
        return bins_dir
    
    # Map reads to contigs for depth file
    command = (
        f"minimap2 -ax sr {contigs} {clean_fastq1} {clean_fastq2} | "
        f"samtools view -bS - | samtools sort -o {bam_file}"
    )
    run_command(command, "Mapping reads for MetaBAT2")
    
    # Create depth file
    depth_file = os.path.join(output_dir, "depth.txt")
    command = f"jgi_summarize_bam_contig_depths --outputDepth {depth_file} {bam_file}"
    run_command(command, "Creating depth file for MetaBAT2")
    
    # Run MetaBAT2
    command = f"metabat2 -i {contigs} -a {depth_file} -o {bins_dir}/bin -v"
    run_command(command, "Binning with MetaBAT2")
    return bins_dir

def check_and_download_kraken2_db(db_path=None):
    """Check for Kraken2 database and download MiniKraken if missing."""
    kraken_db_path = db_path or os.environ.get('KRAKEN2_DB_PATH')
    if not kraken_db_path:
        error_msg = (
            "KRAKEN2_DB_PATH environment variable or --db_path argument not set. "
            "Please set the environment variable with 'export KRAKEN2_DB_PATH=/path/to/database' "
            "or provide --db_path /path/to/database in the command line."
        )
        logging.error(error_msg)
        raise ValueError(error_msg)
    
    db_dir = Path(kraken_db_path)
    # Check for a key database file to confirm database presence
    db_marker = db_dir / "hash.k2d"
    
    if db_marker.exists():
        logging.info("Kraken2 database found")
        return kraken_db_path
    
    logging.info("Kraken2 database not found, downloading MiniKraken...")
    try:
        # Ensure the directory exists
        db_dir.mkdir(parents=True, exist_ok=True)
        # Download and extract MiniKraken2 v2 8GB database
        command = (
            f"wget -P {db_dir} https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz && "
            f"tar -xzf {db_dir}/k2_standard_08gb_20250402.tar.gz -C {db_dir} && "
            f"rm {db_dir}/k2_standard_08gb_20250402.tar.gz"
        )
        run_command(command, "Downloading Kraken2 MiniKraken database")
        logging.info("Kraken2 MiniKraken database downloaded successfully")
    except Exception as e:
        logging.error(f"Failed to download Kraken2 database: {str(e)}")
        raise

def taxonomic_classification(bins_dir, output_dir, resume=False, db_path=None):
    """Perform taxonomic classification using Kraken2 for each bin separately."""
    kraken_output_dir = os.path.join(output_dir, "kraken2_output")
    os.makedirs(kraken_output_dir, exist_ok=True)
    
    # Check and download database if needed
    kraken_db_path = check_and_download_kraken2_db(db_path)
    
    # Get list of bin files
    bin_files = list(Path(bins_dir).glob("bin*.fa"))
    if not bin_files:
        logging.error("No bin files found in {bins_dir}, cannot perform taxonomic classification")
        raise FileNotFoundError(f"No bin files found in {bins_dir}")
    
    kraken_output_files = []
    for bin_file in bin_files:
        bin_id = bin_file.stem
        kraken_output = os.path.join(kraken_output_dir, f"{bin_id}_kraken2_output.txt")
        kraken_report = os.path.join(kraken_output_dir, f"{bin_id}_kraken2_report.txt")
        
        if resume and os.path.exists(kraken_output) and os.path.exists(kraken_report):
            logging.info(f"Kraken2 output files for {bin_id} found, skipping")
        else:
            command = (
                f"kraken2 --db {kraken_db_path} --threads 16 "
                f"--output {kraken_output} --report {kraken_report} --use-names {bin_file}"
            )
            run_command(command, f"Taxonomic classification with Kraken2 for {bin_id}")
        kraken_output_files.append(kraken_output)
    
    if not kraken_output_files:
        logging.error("No Kraken2 output files generated")
        raise RuntimeError("Kraken2 failed to generate output files")
    
    return kraken_output_files

def detect_amr_virulence_abricate(bins_dir, output_dir, resume=False):
    """Detect AMR, virulence, and plasmid genes using ABRicate with CARD, VFDB, PlasmidFinder, and ResFinder databases."""
    abricate_output_dir = os.path.join(output_dir, "abricate")
    databases = {
        "card": os.path.join(abricate_output_dir, "card", "card_summary.tsv"),
        "vfdb": os.path.join(abricate_output_dir, "vfdb", "vfdb_summary.tsv"),
        "plasmidfinder": os.path.join(abricate_output_dir, "plasmidfinder", "plasmidfinder_summary.tsv"),
        "resfinder": os.path.join(abricate_output_dir, "resfinder", "resfinder_summary.tsv")
    }
    
    # Check if all summary files exist for resume
    if resume and all(os.path.exists(summary) for summary in databases.values()):
        logging.info("All ABRicate outputs found, skipping detection")
        return list(databases.values())
    
    # Check if bin files exist
    bin_files = list(Path(bins_dir).glob("bin*.fa"))
    if not bin_files:
        logging.warning("No bin files found for ABRicate, skipping detection")
        return [None] * len(databases)
    
    results = []
    for db_name, summary_file in databases.items():
        db_output_dir = os.path.dirname(summary_file)
        os.makedirs(db_output_dir, exist_ok=True)
        
        # Run ABRicate for each bin
        for bin_file in bin_files:
            bin_id = bin_file.stem
            output_file = os.path.join(db_output_dir, f"{bin_id}_{db_name}.tsv")
            command = (
                f"abricate --db {db_name} --quiet {bin_file} > {output_file}"
            )
            run_command(command, f"{db_name.capitalize()} detection with ABRicate for {bin_id}")
        
        # Summarize results
        command = f"abricate --summary {db_output_dir}/*.tsv > {summary_file}"
        run_command(command, f"Summarizing ABRicate {db_name.capitalize()} results")
        results.append(summary_file)
    
    return results

def map_and_segregate_fastq(contigs, clean_fastq1, clean_fastq2, taxonomy_files, bins_dir, output_dir, resume=False):
    """Map FASTQ to binned contigs and segregate by organism into bin-specific subdirectories."""
    logging.info(f"Using contigs file for mapping: {contigs}")
    # Parse Kraken2 output files to map contigs to taxonomic labels
    taxonomy = {}
    contig_to_taxa = {}
    for taxonomy_file in taxonomy_files:
        bin_id = Path(taxonomy_file).stem.replace('_kraken2_output', '')
        with open(taxonomy_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 3:
                    continue
                classification = fields[0]  # 'C' or 'U' for classified/unclassified
                contig_id = fields[1]
                taxon_name = fields[2]
                if classification == 'C':
                    contig_to_taxa[contig_id] = sanitize_organism_name(taxon_name)
                else:
                    contig_to_taxa[contig_id] = "unclassified"
    
    # Map contigs to bins
    bin_to_contigs = {}
    for bin_file in Path(bins_dir).glob("bin*.fa"):
        bin_id = bin_file.stem
        bin_to_contigs[bin_id] = []
        try:
            with open(bin_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        contig_id = line.strip()[1:].split()[0]  # Take first word to handle spaces
                        bin_to_contigs[bin_id].append(contig_id)
        except Exception as e:
            logging.warning(f"Failed to read bin file {bin_file}: {str(e)}")
            continue
    
    # Assign taxonomy to bins (use most common taxon among contigs)
    for bin_id, contigs_list in bin_to_contigs.items():
        if not contigs_list:
            logging.warning(f"Bin {bin_id} contains no contigs, skipping")
            taxonomy[bin_id] = "unclassified"
            continue
        taxa_counts = Counter(contig_to_taxa.get(contig, "unclassified") for contig in contigs_list)
        most_common_taxa = taxa_counts.most_common(1)
        taxonomy[bin_id] = most_common_taxa[0][0] if most_common_taxa else "unclassified"
        logging.info(f"Bin {bin_id} assigned taxon: {taxonomy[bin_id]}")
    
    # Check if any segregated FASTQ files exist
    all_segregated = True
    for bin_id, organism in taxonomy.items():
        bin_subdir = os.path.join(bins_dir, bin_id)
        output_fastq1 = os.path.join(bin_subdir, f"{organism}_R1.fastq")
        output_fastq2 = os.path.join(bin_subdir, f"{organism}_R2.fastq")
        if not (os.path.exists(output_fastq1) and os.path.exists(output_fastq2)):
            all_segregated = False
            break
    if resume and all_segregated:
        logging.info("All segregated FASTQ files found in bin subdirectories, skipping mapping and segregation")
        return
    
    # Map reads to contigs
    bam_file = os.path.join(output_dir, "mapped.bam")
    command = (
        f"minimap2 -ax sr {contigs} {clean_fastq1} {clean_fastq2} | "
        f"samtools view -bS - | samtools sort -o {bam_file}"
    )
    run_command(command, "Mapping reads to contigs")
    
    # Extract reads per bin
    read_to_bin = {}
    cmd = ["samtools", "view", bam_file]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        fields = line.strip().split('\t')
        if len(fields) < 3:
            continue
        read_id = fields[0]
        contig_id = fields[2]
        if contig_id != "*" and contig_id in contig_to_taxa:
            for bin_id, contigs_list in bin_to_contigs.items():
                if contig_id in contigs_list:
                    read_to_bin[read_id] = bin_id
                    break
    proc.terminate()
    
    # Segregate reads into bin-specific subdirectories
    for bin_id in bin_to_contigs:
        organism = taxonomy.get(bin_id, "unclassified")
        if organism == "unclassified":
            logging.warning(f"Skipping unclassified bin {bin_id}")
            continue
        bin_subdir = os.path.join(bins_dir, bin_id)
        os.makedirs(bin_subdir, exist_ok=True)
        output_fastq1 = os.path.join(bin_subdir, f"{organism}_R1.fastq")
        output_fastq2 = os.path.join(bin_subdir, f"{organism}_R2.fastq")
        
        if resume and os.path.exists(output_fastq1) and os.path.exists(output_fastq2):
            logging.info(f"Segregated FASTQ files for {organism} in bin {bin_id} found, skipping")
            continue
        
        read_ids = {read_id for read_id, bid in read_to_bin.items() if bid == bin_id}
        if not read_ids:
            logging.warning(f"No reads mapped to bin {bin_id} (taxon {organism}), skipping FASTQ generation")
            continue
        
        with open(output_fastq1, 'w') as f1, open(output_fastq2, 'w') as f2:
            for record in SeqIO.parse(clean_fastq1, "fastq"):
                if record.id in read_ids:
                    SeqIO.write(record, f1, "fastq")
            for record in SeqIO.parse(clean_fastq2, "fastq"):
                if record.id in read_ids:
                    SeqIO.write(record, f2, "fastq")
        
        # Verify FASTQ files
        if os.path.getsize(output_fastq1) == 0 or os.path.getsize(output_fastq2) == 0:
            logging.warning(f"Empty FASTQ files generated for bin {bin_id} (taxon {organism}), removing")
            os.remove(output_fastq1)
            os.remove(output_fastq2)
        else:
            logging.info(f"Generated FASTQ files for bin {bin_id} (taxon {organism}): {output_fastq1}, {output_fastq2}")

def main(input_fastq1, input_fastq2, output_dir, resume=False, db_path=None):
    """Main pipeline function."""
    try:
        os.makedirs(output_dir, exist_ok=True)
        logging.info("Starting metagenomics pipeline")
        
        # Quality control
        clean_fastq1, clean_fastq2 = quality_control(input_fastq1, input_fastq2, output_dir, resume)
        
        # De novo assembly
        contigs = denovo_assembly(clean_fastq1, clean_fastq2, output_dir, resume)
        
        # Binning
        bins_dir = bin_contigs(contigs, clean_fastq1, clean_fastq2, output_dir, resume)
        
        # Taxonomic classification
        taxonomy_files = taxonomic_classification(bins_dir, output_dir, resume, db_path)
        
        # AMR, virulence, and plasmid detection with ABRicate
        card_summary, vfdb_summary, plasmidfinder_summary, resfinder_summary = detect_amr_virulence_abricate(bins_dir, output_dir, resume)
        
        # Map and segregate FASTQ
        map_and_segregate_fastq(contigs, clean_fastq1, clean_fastq2, taxonomy_files, bins_dir, output_dir, resume)
        
        logging.info("Pipeline completed successfully")
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Metagenomics pipeline with resume capability")
    parser.add_argument("--fastq1", required=True, help="Path to forward FASTQ file")
    parser.add_argument("--fastq2", required=True, help="Path to reverse FASTQ file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--resume", action="store_true", help="Resume pipeline from last completed step")
    parser.add_argument("--db_path", help="Path to Kraken2 database (overrides KRAKEN2_DB_PATH)")
    args = parser.parse_args()
    
    main(args.fastq1, args.fastq2, args.output, args.resume, args.db_path)