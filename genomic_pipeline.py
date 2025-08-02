#!/usr/bin/env python3
import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

def run_command(cmd_list, error_msg, output_file=None, input_file=None, pipe_to=None, binary=False):
    """Run a command, optionally with input file, output file, or piped command, handling binary or text output."""
    try:
        if pipe_to:
            p1 = subprocess.run(
                cmd_list,
                check=True,
                text=not binary,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                input=open(input_file, 'rb').read() if input_file and binary else open(input_file, 'r').read() if input_file else None
            )
            p2 = subprocess.run(
                pipe_to,
                check=True,
                text=not binary,
                stdout=open(output_file, 'wb') if output_file and binary else open(output_file, 'w') if output_file else subprocess.PIPE,
                stderr=subprocess.PIPE,
                input=p1.stdout
            )
            if p2.stderr:
                logging.debug(f"Second command stderr: {p2.stderr.decode('utf-8') if binary else p2.stderr}")
        else:
            subprocess.run(
                cmd_list,
                check=True,
                text=not binary,
                stdout=open(output_file, 'wb') if output_file and binary else open(output_file, 'w') if output_file else None,
                stderr=subprocess.PIPE,
                input=open(input_file, 'rb').read() if input_file and binary else open(input_file, 'r').read() if input_file else None
            )
    except subprocess.CalledProcessError as e:
        stderr = e.stderr.decode('utf-8') if binary else e.stderr
        logging.error(f"{error_msg}\n{stderr}")
        sys.exit(1)

def check_file(file_path, file_desc):
    if not os.path.isfile(file_path):
        logging.error(f"{file_desc} ({file_path}) not found.")
        sys.exit(1)

def validate_fastq(file_path):
    """Check if a file is in valid FASTQ format."""
    try:
        with open(file_path, 'r') as f:
            lines = [f.readline().strip() for _ in range(4)]
            if len(lines) < 4:
                logging.error(f"File {file_path} is too short to be a valid FASTQ file.")
                sys.exit(1)
            if not (lines[0].startswith('@') and lines[2].startswith('+') and len(lines[1]) == len(lines[3])):
                logging.error(f"File {file_path} does not conform to FASTQ format.")
                sys.exit(1)
        logging.info(f"File {file_path} passed FASTQ format validation.")
    except Exception as e:
        logging.error(f"Error validating FASTQ file {file_path}: {e}")
        sys.exit(1)

def validate_bam(bam_file):
    """Check if a BAM file is valid and indexed."""
    try:
        check_file(bam_file, "BAM file")
        if os.path.getsize(bam_file) == 0:
            logging.error(f"BAM file {bam_file} is empty.")
            sys.exit(1)
        subprocess.run(["samtools", "quickcheck", bam_file], check=True, stderr=subprocess.PIPE)
        bam_index = bam_file + ".bai"
        if not os.path.exists(bam_index):
            logging.info(f"Indexing BAM file {bam_file}")
            run_command(["samtools", "index", bam_file], f"Failed to index BAM file {bam_file}")
        logging.info(f"BAM file {bam_file} passed validation.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Invalid BAM file {bam_file}: {e.stderr.decode('utf-8')}")
        sys.exit(1)

def index_reference_fasta(fasta_file):
    """Index a FASTA file for bcftools if not already indexed."""
    fai_file = fasta_file + ".fai"
    if not os.path.exists(fai_file):
        logging.info(f"Indexing reference FASTA file {fasta_file}")
        run_command(["samtools", "faidx", fasta_file], f"Failed to index reference FASTA {fasta_file}")
    else:
        logging.info(f"Reference FASTA index {fai_file} already exists.")

def validate_genbank(genbank_file):
    """Validate GenBank file."""
    try:
        check_file(genbank_file, "GenBank file")
        with open(genbank_file, 'r') as f:
            content = f.read()
            if not content.startswith("LOCUS"):
                logging.error(f"Invalid GenBank file {genbank_file}: Does not start with 'LOCUS'")
                sys.exit(1)
            if "FEATURES" not in content:
                logging.error(f"GenBank file {genbank_file}: Lacks 'FEATURES' section. Annotations may be missing.")
                sys.exit(1)
        logging.info(f"GenBank file {genbank_file} passed validation.")
    except Exception as e:
        logging.error(f"Error validating GenBank file {genbank_file}: {e}")
        sys.exit(1)

def validate_fasta_genbank_match(genbank_file, reference_fasta):
    """Validate that sequence IDs in FASTA and GenBank files match."""
    try:
        with open(reference_fasta, "r") as f:
            fasta_id = f.readline().strip().lstrip(">").split()[0]
        with open(genbank_file, "r") as f:
            for line in f:
                if line.startswith("LOCUS"):
                    gbk_id = line.split()[1]
                    break
            else:
                logging.error("No LOCUS line found in GenBank file.")
                sys.exit(1)
        if fasta_id != gbk_id:
            logging.error(f"Sequence ID mismatch: FASTA ID '{fasta_id}' does not match GenBank ID '{gbk_id}'.")
            sys.exit(1)
        logging.info(f"Sequence IDs match: {fasta_id}")
    except Exception as e:
        logging.error(f"FASTA-GenBank validation failed: {e}")
        sys.exit(1)

def check_download_eggnog_database():
    """Check for eggNOG DIAMOND database and download using download_eggnog_data.py if missing."""
    # Determine eggNOG data directory: use EGGNOG_DATA_DIR or default to eggnog-mapper's data/ directory
    eggnog_data_dir = os.environ.get("EGGNOG_DATA_DIR")
    if not eggnog_data_dir:
        # Default to eggnog-mapper's data/ directory in Conda environment
        eggnog_data_dir = os.path.join(os.environ.get("CONDA_PREFIX", ""), "lib", "python3.10", "site-packages", "data")
        if not os.path.exists(eggnog_data_dir):
            # If data/ directory doesn't exist, try eggnog-mapper's root data/ directory
            eggnog_mapper_dir = os.path.join(os.environ.get("CONDA_PREFIX", ""), "lib", "python3.10", "site-packages", "eggnogmapper")
            eggnog_data_dir = os.path.join(eggnog_mapper_dir, "data")
    
    eggnog_db = os.path.join(eggnog_data_dir, "eggnog_proteins.dmnd")
    
    if os.path.exists(eggnog_db) and os.path.getsize(eggnog_db) > 0:
        logging.info(f"eggNOG DIAMOND database found at {eggnog_db}")
        return True
    
    logging.info(f"eggNOG DIAMOND database not found at {eggnog_db}. Attempting to download using download_eggnog_data.py...")
    os.makedirs(eggnog_data_dir, exist_ok=True)
    
    download_script = os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin", "download_eggnog_data.py")
    if not os.path.exists(download_script):
        logging.error(f"download_eggnog_data.py not found at {download_script}. Please ensure eggNOG-mapper is installed correctly.")
        logging.error("To download manually:")
        logging.error("1. Activate Conda environment: 'source activate genomic_pipeline'")
        logging.error(f"2. Run: 'python {download_script}'")
        logging.error(f"3. Or download 'eggnog_proteins.dmnd.gz' from http://eggnog5.embl.de/#/app/downloads")
        logging.error(f"4. Decompress and place it at {eggnog_db}")
        return False
    
    cmd = [sys.executable, download_script, "--data_dir", eggnog_data_dir]
    try:
        logging.info(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
        if os.path.exists(eggnog_db) and os.path.getsize(eggnog_db) > 0:
            logging.info(f"Successfully downloaded eggNOG DIAMOND database to {eggnog_db}")
            return True
        else:
            logging.error(f"Download completed but {eggnog_db} is missing or empty")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to download eggNOG DIAMOND database: {e.stderr}")
    
    logging.error("Manual download instructions:")
    logging.error("1. Activate Conda environment: 'source activate genomic_pipeline'")
    logging.error(f"2. Run: 'python {download_script} --data_dir {eggnog_data_dir}'")
    logging.error(f"3. Or download 'eggnog_proteins.dmnd.gz' from http://eggnog5.embl.de/#/app/downloads")
    logging.error(f"4. Decompress and place it at {eggnog_db}")
    return False

def pair_fastq(fastq1, fastq2, output_fastq1, output_fastq2):
    """Filter paired-end FASTQ files to keep only matched reads."""
    logging.info(f"Pairing FASTQ files {fastq1} and {fastq2}")
    paired_records = []
    try:
        with open(fastq1, 'r') as f1, open(fastq2, 'r') as f2:
            r1_records = {}
            r2_records = {}
            while True:
                r1_lines = [f1.readline().strip() for _ in range(4)]
                if not r1_lines[0]:
                    break
                if not (r1_lines[0].startswith('@') and r1_lines[2].startswith('+') and len(r1_lines[1]) == len(r1_lines[3])):
                    logging.error(f"Invalid FASTQ format in {fastq1} at record starting with {r1_lines[0]}")
                    sys.exit(1)
                seq_id = r1_lines[0][1:].split()[0].rstrip('/1')
                r1_records[seq_id] = r1_lines
            while True:
                r2_lines = [f2.readline().strip() for _ in range(4)]
                if not r2_lines[0]:
                    break
                if not (r2_lines[0].startswith('@') and r2_lines[2].startswith('+') and len(r2_lines[1]) == len(r2_lines[3])):
                    logging.error(f"Invalid FASTQ format in {fastq2} at record starting with {r2_lines[0]}")
                    sys.exit(1)
                seq_id = r2_lines[0][1:].split()[0].rstrip('/2')
                r2_records[seq_id] = r2_lines
            common_ids = set(r1_records.keys()) & set(r2_records.keys())
            if not common_ids:
                logging.error(f"No matching sequence IDs found between {fastq1} and {fastq2}")
                sys.exit(1)
            for seq_id in common_ids:
                paired_records.append((seq_id, r1_records[seq_id], r2_records[seq_id]))
            logging.info(f"Found {len(common_ids)} paired reads")
        
        with open(output_fastq1, 'w') as f1, open(output_fastq2, 'w') as f2:
            for _, r1_lines, r2_lines in paired_records:
                for line in r1_lines:
                    f1.write(line + '\n')
                for line in r2_lines:
                    f2.write(line + '\n')
        
        logging.info(f"Paired FASTQ files written to {output_fastq1} and {output_fastq2}")
        validate_fastq(output_fastq1)
        validate_fastq(output_fastq2)
        return len(common_ids)
    except Exception as e:
        logging.error(f"Failed to pair FASTQ files {fastq1} and {fastq2}: {e}")
        sys.exit(1)

def sort_fastq_pair(fastq1, fastq2, output_fastq1, output_fastq2):
    """Sort paired-end FASTQ files by sequence ID, ensuring paired reads are kept in sync."""
    logging.info(f"Sorting paired FASTQ files {fastq1} and {fastq2} by sequence ID")
    records = []
    try:
        with open(fastq1, 'r') as f1, open(fastq2, 'r') as f2:
            while True:
                r1_lines = [f1.readline().strip() for _ in range(4)]
                r2_lines = [f2.readline().strip() for _ in range(4)]
                if not r1_lines[0] or not r2_lines[0]:
                    if r1_lines[0] or r2_lines[0]:
                        logging.error(f"Mismatch in number of reads between {fastq1} and {fastq2}")
                        sys.exit(1)
                    break
                for lines, file_path in [(r1_lines, fastq1), (r2_lines, fastq2)]:
                    if not (lines[0].startswith('@') and lines[2].startswith('+') and len(lines[1]) == len(lines[3])):
                        logging.error(f"Invalid FASTQ format in {file_path} at record starting with {lines[0]}")
                        sys.exit(1)
                seq_id1 = r1_lines[0][1:].split()[0].rstrip('/1')
                seq_id2 = r2_lines[0][1:].split()[0].rstrip('/2')
                if seq_id1 != seq_id2:
                    logging.error(f"Mismatched sequence IDs: {seq_id1} in {fastq1} vs {seq_id2} in {fastq2}")
                    sys.exit(1)
                records.append((seq_id1, r1_lines, r2_lines))
        
        records.sort(key=lambda x: x[0])
        
        with open(output_fastq1, 'w') as f1, open(output_fastq2, 'w') as f2:
            for _, r1_lines, r2_lines in records:
                for line in r1_lines:
                    f1.write(line + '\n')
                for line in r2_lines:
                    f2.write(line + '\n')
        
        logging.info(f"Sorted FASTQ files written to {output_fastq1} and {output_fastq2}")
        validate_fastq(output_fastq1)
        validate_fastq(output_fastq2)
    except Exception as e:
        logging.error(f"Failed to sort paired FASTQ files {fastq1} and {fastq2}: {e}")
        sys.exit(1)

def create_snpeff_config(output_dir, snpeff_db, reference_fasta):
    """Create SnpEff configuration file."""
    snpeff_data_dir = os.path.abspath(os.path.join(output_dir, "snpeff_data"))
    snpeff_config = os.path.join(output_dir, "snpEff.config")
    config_content = f"""# SnpEff configuration file
data.dir = {snpeff_data_dir}

# Database entry
{snpeff_db}.genome = {snpeff_db}
{snpeff_db}.reference = {os.path.abspath(reference_fasta)}
"""
    if not os.path.exists(snpeff_config):
        os.makedirs(output_dir, exist_ok=True)
        with open(snpeff_config, "w") as f:
            f.write(config_content)
        logging.info(f"Created SnpEff configuration file at {snpeff_config}")
    else:
        logging.info(f"SnpEff configuration file already exists at {snpeff_config}")
    return snpeff_config, snpeff_data_dir

def snpeff_entry_exists(config_path, db_name):
    """Check if database entry exists in SnpEff config."""
    try:
        with open(config_path, "r") as f:
            for line in f:
                if line.strip().startswith(f"{db_name}.genome"):
                    return True
    except Exception as e:
        logging.error(f"Error reading SnpEff config file {config_path}: {e}")
        return False
    return False

def check_create_snpeff_db(snpeff_config, snpeff_db, genbank_file, reference_fasta, snpeff_data_dir):
    """Create or validate SnpEff database."""
    db_path = os.path.join(snpeff_data_dir, snpeff_db)
    predictor_file = os.path.join(db_path, "snpEffectPredictor.bin")
    
    validate_genbank(genbank_file)
    validate_fasta_genbank_match(genbank_file, reference_fasta)
    
    if os.path.exists(predictor_file) and snpeff_entry_exists(snpeff_config, snpeff_db):
        logging.info(f"SnpEff database ({snpeff_db}) exists and appears valid.")
        return
    
    logging.info(f"Creating or recreating SnpEff database {snpeff_db}")
    os.makedirs(db_path, exist_ok=True)
    genbank_dest = os.path.join(db_path, "genes.gbk")
    fasta_dest = os.path.join(db_path, "sequences.fa")
    run_command(["cp", genbank_file, genbank_dest], f"Failed to copy GenBank file to {genbank_dest}")
    run_command(["cp", reference_fasta, fasta_dest], f"Failed to copy FASTA file to {fasta_dest}")
    
    if not snpeff_entry_exists(snpeff_config, snpeff_db):
        with open(snpeff_config, "a") as f:
            f.write(f"\n{snpeff_db}.genome : {snpeff_db}\n{snpeff_db}.reference = {os.path.abspath(reference_fasta)}\n")
    
    build_cmd = ["snpEff", "build", "-genbank", "-v", snpeff_db, "-c", snpeff_config, "-dataDir", snpeff_data_dir]
    run_command(build_cmd, f"Failed to create SnpEff database ({snpeff_db})")
    if not os.path.exists(predictor_file):
        logging.error(f"SnpEff database build failed: {predictor_file} not found.")
        sys.exit(1)
    logging.info(f"SnpEff database {snpeff_db} created successfully.")

def mark_step_complete(step, checkpoint_file):
    with open(checkpoint_file, "a") as f:
        f.write(f"{step}\n")

def is_step_complete(step, checkpoint_file):
    if not os.path.exists(checkpoint_file):
        return False
    with open(checkpoint_file, "r") as f:
        completed = f.read().splitlines()
    return step in completed

def main():
    parser = argparse.ArgumentParser(description="Variant calling and annotation pipeline for paired-end reads")
    parser.add_argument("-s", "--sample", required=True, help="Sample name")
    parser.add_argument("-1", "--fastq1", required=True, help="Path to the first paired-end FASTQ file (R1)")
    parser.add_argument("-2", "--fastq2", required=True, help="Path to the second paired-end FASTQ file (R2)")
    parser.add_argument("-r", "--ref_dir", required=True, help="Directory containing reference.fasta and reference.gbk")
    parser.add_argument("-d", "--snpeff_db", required=True, help="SnpEff database name")
    parser.add_argument("-o", "--output_dir", default="output", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of CPU threads")
    parser.add_argument("-k", "--organism", default="eco", help="Organism code for KEGG")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument("--resume", action="store_true", help="Resume from completed steps")
    parser.add_argument("--skip-emapper", action="store_true", help="Skip eggNOG-mapper step if database is unavailable")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO)

    checkpoint_file = os.path.join(args.output_dir, f"{args.sample}.checkpoint")

    if not args.resume and os.path.exists(checkpoint_file):
        logging.info(f"Clearing checkpoint file {checkpoint_file} for fresh run")
        os.remove(checkpoint_file)

    reference = os.path.join(args.ref_dir, "reference.fasta")
    genbank_file = os.path.join(args.ref_dir, "reference.gbk")
    for f, desc in [(args.fastq1, "FASTQ R1"), (args.fastq2, "FASTQ R2"), (reference, "Reference FASTA"), (genbank_file, "GenBank")]:
        check_file(f, desc)

    index_reference_fasta(reference)

    os.makedirs(args.output_dir, exist_ok=True)
    snpeff_config, snpeff_data_dir = create_snpeff_config(args.output_dir, args.snpeff_db, reference)

    sample = args.sample
    paired_fastq1 = os.path.join(args.output_dir, f"{sample}_R1.paired.fastq")
    paired_fastq2 = os.path.join(args.output_dir, f"{sample}_R2.paired.fastq")
    sorted_fastq1 = os.path.join(args.output_dir, f"{sample}_R1.sorted.fastq")
    sorted_fastq2 = os.path.join(args.output_dir, f"{sample}_R2.sorted.fastq")
    sam = os.path.join(args.output_dir, f"{sample}.sam")
    bam = os.path.join(args.output_dir, f"{sample}.sorted.bam")
    vcf = os.path.join(args.output_dir, f"{sample}.vcf.gz")
    consensus = os.path.join(args.output_dir, f"{sample}.consensus.fasta")
    filtered_vcf = os.path.join(args.output_dir, f"{sample}.filtered.vcf.gz")
    annotated_vcf = os.path.join(args.output_dir, f"{sample}.annotated.vcf")
    prokka_dir = os.path.join(args.output_dir, "prokka")
    emapper_out = os.path.join(args.output_dir, sample)
    annot_file = os.path.join(args.output_dir, f"{sample}.annotations.txt")

    check_create_snpeff_db(snpeff_config, args.snpeff_db, genbank_file, reference, snpeff_data_dir)
    if not is_step_complete("create_snpeff_db", checkpoint_file):
        mark_step_complete("create_snpeff_db", checkpoint_file)

    has_eggnog_db = check_download_eggnog_database()
    if not has_eggnog_db and not args.skip_emapper:
        logging.error("eggNOG DIAMOND database is required for the pipeline. Use --skip-emapper to proceed without eggNOG-mapper.")
        sys.exit(1)
    elif not has_eggnog_db:
        logging.warning("Skipping eggNOG-mapper step due to missing database and --skip-emapper flag.")

    pipeline_steps = [
        ("pair_fastq", None, "Pair FASTQ files to keep matched reads", (paired_fastq1, paired_fastq2), None, None, False),
        ("sort_fastq", None, "Sort paired FASTQ files by sequence ID", (sorted_fastq1, sorted_fastq2), None, None, False),
        ("index_reference", ["bowtie2-build", "--threads", str(args.threads), reference, os.path.join(args.ref_dir, "reference")], "Index reference genome with bowtie2", None, None, None, False),
        ("mapping", ["bowtie2", "-x", os.path.join(args.ref_dir, "reference"), "-1", sorted_fastq1, "-2", sorted_fastq2, "-p", str(args.threads)], f"Map paired-end reads > {sam}", sam, None, None, False),
        ("sam_to_bam_sort", ["samtools", "view", "-bS"], f"Convert SAM to sorted BAM > {bam}", bam, sam, ["samtools", "sort", "-o", bam], True),
        ("index_bam", ["samtools", "index", bam], "Index BAM", None, None, None, False),
        ("mpileup_call_variants", ["bcftools", "mpileup", "-Ou", "-f", reference, bam], f"Run mpileup and call variants > {vcf}", vcf, None, ["bcftools", "call", "-mv", "-Oz", "-o", vcf], True),
        ("index_vcf", ["bcftools", "index", vcf], "Index VCF", None, None, None, False),
        ("make_consensus", ["bcftools", "consensus", "-f", reference, vcf, "-o", consensus], "Create consensus", consensus, None, None, False),
        ("filter_vcf", ["bcftools", "filter", "-i", "QUAL>20 && DP>10", vcf, "-Oz", "-o", filtered_vcf], "Filter VCF", filtered_vcf, None, None, True),
        ("annotate_snpeff", ["snpEff", "-v", args.snpeff_db, "-c", snpeff_config, "-dataDir", snpeff_data_dir, filtered_vcf], f"Annotate VCF > {annotated_vcf}", annotated_vcf, None, None, False),
        ("run_prokka", ["prokka", "--outdir", prokka_dir, "--prefix", sample, "--cpus", str(args.threads), consensus], "Run Prokka", None, None, None, False),
    ]

    if has_eggnog_db and not args.skip_emapper:
        pipeline_steps.append(
            ("run_emapper", ["emapper.py", "-i", os.path.join(prokka_dir, f"{sample}.faa"), "-o", emapper_out, "--cpu", str(args.threads), "-m", "diamond"], "Run eggNOG-mapper", None, None, None, False)
        )

    for step, cmd, desc, output_file, input_file, pipe_to, binary in pipeline_steps:
        if args.resume and is_step_complete(step, checkpoint_file):
            logging.info(f"Skipping completed step: {step}")
            continue
        if args.dry_run:
            if step in ["pair_fastq", "sort_fastq"]:
                logging.info(f"[Dry-run] {desc}: Python-based {'pairing' if step == 'pair_fastq' else 'paired FASTQ sorting'}")
            elif pipe_to:
                logging.info(f"[Dry-run] {desc}: {' '.join(cmd)} < {input_file if input_file else ''} | {' '.join(pipe_to)}")
            else:
                logging.info(f"[Dry-run] {desc}: {' '.join(cmd)}")
        else:
            logging.info(f"Running: {desc}")
            if step == "pair_fastq":
                pair_fastq(args.fastq1, args.fastq2, paired_fastq1, paired_fastq2)
            elif step == "sort_fastq":
                sort_fastq_pair(paired_fastq1, paired_fastq2, sorted_fastq1, sorted_fastq2)
            elif step == "index_bam":
                validate_bam(bam)
                run_command(cmd, f"Failed at: {desc}", output_file, input_file, pipe_to, binary)
            elif step == "mpileup_call_variants":
                validate_bam(bam)
                run_command(cmd, f"Failed at: {desc}", output_file, input_file, pipe_to, binary)
            else:
                run_command(cmd, f"Failed at: {desc}", output_file, input_file, pipe_to, binary)
        mark_step_complete(step, checkpoint_file)

    step = "extract_annotations"
    if (has_eggnog_db and not args.skip_emapper) and (not args.resume or not is_step_complete(step, checkpoint_file)):
        try:
            with open(f"{emapper_out}.emapper.annotations", "r") as f:
                # Find the last line starting with '#'
                header_line = None
                for line in f:
                    if line.startswith("#"):
                        header_line = line.strip()
                    else:
                        break  # Stop at the first non-comment line to optimize
                f.seek(0)  # Reset file pointer to start
                # Define the column indices to extract
                selected_indices = [0, 6, 7, 8, 9, 11, 12, 13, 14, 15]
                # Map indices to desired column names for compatibility with assess_pathway_completeness.py
                header_mapping = {
                    0: "Gene",
                    6: "Description",
                    7: "Preferred_name",
                    8: "GO",
                    9: "EC",
                    11: "KEGG",
                    12: "KEGG_Pathway",
                    13: "KEGG_Module",
                    14: "KEGG_Reaction",
                    15: "KEGG_rclass"
                }
                with open(annot_file, "w") as out:
                    # Write header with renamed columns
                    if header_line:
                        header_fields = header_line.split("\t")
                        selected_headers = [header_mapping.get(i, header_fields[i].lstrip("#")) for i in selected_indices]
                        out.write("\t".join(selected_headers) + "\n")
                    else:
                        logging.warning("No header line found in eggNOG-mapper annotations file. Writing data without header.")
                    # Write data rows
                    for line in f:
                        if line.startswith("#"): 
                            continue
                        fields = line.strip().split("\t")
                        if len(fields) >= 12 and (fields[8] or fields[11]):
                            selected_fields = [fields[i] for i in selected_indices if i < len(fields)]
                            out.write("\t".join(selected_fields) + "\n")
            mark_step_complete(step, checkpoint_file)
        except FileNotFoundError:
            logging.warning(f"eggNOG-mapper annotations file {emapper_out}.emapper.annotations not found. Skipping annotation extraction.")
    else:
        logging.info(f"Skipping completed or disabled step: {step}")

    if not args.keep_temp:
        for temp_file in [paired_fastq1, paired_fastq2, sorted_fastq1, sorted_fastq2, sam]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                logging.info(f"Removed temporary file: {temp_file}")

    logging.info("Pipeline completed successfully.")

if __name__ == "__main__":
    main()