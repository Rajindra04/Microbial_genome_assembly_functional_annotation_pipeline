# Metagenomics and Genomic Analysis Pipeline

This pipeline performs metagenomic and genomic analysis on paired-end FASTQ files, processing them through quality control, de novo assembly, binning, taxonomic classification, and functional annotation. It is designed for metagenomic datasets, producing bin-specific FASTQ files and annotations, followed by genomic analysis and visualization for each bin.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Setup Instructions](#setup-instructions)
- [Usage](#usage)
  - [Workflow](#workflow)
  - [Running the Metagenomics Pipeline](#running-the-metagenomics-pipeline)
  - [Running the Genomic Pipeline](#running-the-genomic-pipeline)
  - [Running the Visualization Pipeline](#running-the-visualization-pipeline)
  - [Example Commands](#example-commands)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Overview

The pipeline consists of three main components:
1. **Metagenomics Pipeline** (`metagenomics_pipeline.py`): Processes paired-end FASTQ files through quality control (fastp), de novo assembly (SPAdes), binning (MetaBAT2), taxonomic classification (Kraken2), and AMR/virulence/plasmid detection (ABRicate). It segregates reads into bin-specific FASTQ files based on taxonomic assignments.
2. **Genomic Pipeline** (`genomic_pipeline.py`): Performs variant calling, annotation (Prokka, eggNOG-mapper, SnpEff), and functional analysis on bin-specific FASTQ files, requiring organism-specific reference genomes.
3. **Visualization Pipeline** (`visualize.py`): Generates visualizations (e.g., COG distribution, pathway heatmaps) from genomic annotations.

The pipeline is orchestrated by `main.py`, which supports three modes: `metagenomics`, `genomic`, and `visualize`. The `metagenomics` mode must be run first, followed by `genomic` and `visualize` for each bin, with appropriate reference genomes specified for the `genomic` mode.

## Installation

### Prerequisites
- **Operating System**: Linux or macOS (Windows may require WSL2).
- **Python**: Version 3.8 or higher.
- **Dependencies**:
  - **Bioinformatics Tools**:
    - `fastp` (v0.23.2+): For quality control.
    - `spades.py` (SPAdes v3.15.5+): For de novo assembly.
    - `minimap2` (v2.24+): For read mapping.
    - `samtools` (v1.15+): For BAM file processing.
    - `metabat2` (v2.15+): For binning.
    - `jgi_summarize_bam_contig_depths` (MetaBAT2 suite): For depth calculation.
    - `kraken2` (v2.1.2+): For taxonomic classification.
    - `abricate` (v1.0.1+): For AMR/virulence/plasmid detection.
    - `prokka` (v1.14.6+): For genome annotation.
    - `diamond` (v2.0.15+): For eggNOG-mapper.
    - `snpEff` (v5.0+): For variant annotation.
  - **Python Packages**:
    - `biopython` (v1.79+): For FASTQ parsing.
- **Databases**:
  - **Kraken2 Database**: A Kraken2 database (e.g., MiniKraken, Standard, or PlusPFP) for taxonomic classification.
  - **eggNOG Database**: Required for functional annotation in `eggnog-mapper`.
  - **Reference Genomes**: Organism-specific `reference.fasta` and `reference.gbk` files for each bin.
  - **SnpEff Database**: Organism-specific SnpEff database for variant annotation.

### Setup Instructions
## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/metagenomics_analysis.git
   cd metagenomics_analysis
   ```

2. **Set Up the Conda Environment**:
   ```bash
   chmod +x build.sh
   ./build.sh
   ```

   This script creates a Conda environment named `metagenomics_pipeline` and installs all dependencies listed in `environment.yml`.

3. **Activate the Environment**:
   ```bash
   conda activate metagenomics_pipeline
   ```

4. **Set Up Kraken2 Database**:
   - Specify the path to an existing Kraken2 database via the `KRAKEN2_DB_PATH` environment variable:
     ```bash
     export KRAKEN2_DB_PATH=/path/to/kraken2/database
     ```
   - Alternatively, provide the `--db_path` argument when running the pipeline:
     ```bash
     python main.py --pipeline all --db_path /path/to/kraken2/database ...
     ```
   - If no database is found, the pipeline will attempt to download the MiniKraken2 database (requires internet access).

5. **Set Up eggNOG Database**:
   - Download the eggNOG 5.0 database for `eggnog-mapper`:
     ```bash
     export EGGNOG_DATA_DIR=/path/to/eggnog_db
     mkdir -p $EGGNOG_DATA_DIR
     wget -P $EGGNOG_DATA_DIR http://eggnog5.embl.de/download/eggnog_5.0/eggnog.db.gz
     gunzip $EGGNOG_DATA_DIR/eggnog.db.gz
     ```
   - Alternatively, use `eggnog-mapper`’s built-in download command:
     ```bash
     download_eggnog_data.py -y -P $EGGNOG_DATA_DIR
     ```
   - Ensure `EGGNOG_DATA_DIR` is set or configure `eggnog-mapper` to use the database path:
     ```bash
     export EGGNOG_DATA_DIR=/path/to/eggnog_db
     ```

6. **Prepare Reference Genomes**:
   - Create a directory for each organism (e.g., `ref_dir_Lactiplantibacillus_plantarum/`) containing `reference.fasta` and `reference.gbk`.
   - Example:
     ```
     ref_dir_Lactiplantibacillus_plantarum/
     ├── reference.fasta
     ├── reference.gbk
     ref_dir_Salmonella_enterica/
     ├── reference.fasta
     ├── reference.gbk
     ```

7. **Set Up SnpEff Database**:
   - Configure organism-specific SnpEff databases (refer to SnpEff documentation).

## Usage

### Workflow
1. **Run Metagenomics Pipeline**:
   - Processes input FASTQ files to generate bin-specific FASTQ files and Kraken2 taxonomic outputs.
   - Outputs are stored in `output_dir/bins/` and `output_dir/kraken2_output/`.
   - Check `pipeline.log` for bin-specific taxonomic assignments to identify appropriate reference genomes.

2. **Run Genomic Pipeline**:
   - Run for each bin separately, specifying a bin ID and organism-specific reference directory.
   - Outputs are stored in `output_dir/genomic/binID_organism/`.

3. **Run Visualization Pipeline**:
   - Run for a specific bin or all bins to generate visualizations from genomic annotations.
   - Outputs are stored in `output_dir/visualizations/binID_organism/`.

### Running the Metagenomics Pipeline
```bash
python main.py --pipeline metagenomics \
    --fastq1 input_R1.fastq \
    --fastq2 input_R2.fastq \
    --output output_dir \
    --db_path /path/to/kraken2/database \
    --threads 8 \
    --resume
```
- **Required Arguments**:
  - `--fastq1`, `--fastq2`: Paired-end FASTQ files.
  - `--output`: Output directory.
  - `--db_path`: Path to Kraken2 database (or set `KRAKEN2_DB_PATH` environment variable).
- **Optional Arguments**:
  - `--threads`: Number of CPU threads (default: 8).
  - `--resume`: Skip completed steps if outputs exist.

### Running the Genomic Pipeline
Run for each bin with the appropriate reference directory:
```bash
python main.py --pipeline genomic \
    --bin_id bin1 \
    --ref_dir ref_dir_Lactiplantibacillus_plantarum \
    --snpeff_db lactiplantibacillus_plantarum_db \
    --output output_dir \
    --threads 8 \
    --resume
```
- **Required Arguments**:
  - `--bin_id`: Bin ID (e.g., `bin1`).
  - `--ref_dir`: Directory containing `reference.fasta` and `reference.gbk`.
  - `--snpeff_db`: SnpEff database name.
  - `--output`: Output directory (same as metagenomics pipeline).
- **Optional Arguments**:
  - `--threads`, `--resume`, `--keep-temp`, `--dry-run`, `--verbose`, `--skip-emapper`, `--organism`.

### Running the Visualization Pipeline
Run for a specific bin or all bins:
```bash
python main.py --pipeline visualize \
    --bin_id bin1 \
    --output output_dir
```
- **Required Arguments**:
  - `--output`: Output directory.
- **Optional Arguments**:
  - `--bin_id`: Bin ID to process a single bin (omit to process all genomic outputs).

### Example Commands
1. **Metagenomics Pipeline**:
   ```bash
   python main.py --pipeline metagenomics \
       --fastq1 input_R1.fastq \
       --fastq2 input_R2.fastq \
       --output output_dir \
       --db_path /path/to/kraken2_db \
       --threads 8 \
       --resume
   ```

2. **Genomic Pipeline for bin1**:
   ```bash
   python main.py --pipeline genomic \
       --bin_id bin1 \
       --ref_dir ref_dir_Lactiplantibacillus_plantarum \
       --snpeff_db lactiplantibacillus_plantarum_db \
       --output output_dir \
       --threads 8 \
       --resume
   ```

3. **Genomic Pipeline for bin2**:
   ```bash
   python main.py --pipeline genomic \
       --bin_id bin2 \
       --ref_dir ref_dir_Salmonella_enterica \
       --snpeff_db salmonella_enterica_db \
       --output output_dir \
       --threads 8 \
       --resume
   ```

4. **Visualization for bin1**:
   ```bash
   python main.py --pipeline visualize \
       --bin_id bin1 \
       --output output_dir
   ```

5. **Visualization for all bins**:
   ```bash
   python main.py --pipeline visualize \
       --output output_dir
   ```

## Directory Structure
```
output_dir/
├── bins/
│   ├── bin1/
│   │   ├── bin1.fa
│   │   ├── Lactiplantibacillus_plantarum_R1.fastq
│   │   ├── Lactiplantibacillus_plantarum_R2.fastq
│   ├── bin2/
│   │   ├── bin2.fa
│   │   ├── Salmonella_enterica_R1.fastq
│   │   ├── Salmonella_enterica_R2.fastq
│   └── ...
├── kraken2_output/
│   ├── bin1_kraken2_output.txt
│   ├── bin1_kraken2_report.txt
│   ├── bin2_kraken2_output.txt
│   ├── bin2_kraken2_report.txt
│   └── ...
├── spades_output/
│   ├── contigs.fasta
├── genomic/
│   ├── bin1_Lactiplantibacillus_plantarum/
│   │   ├── bin1_Lactiplantibacillus_plantarum.annotations.txt
│   │   ├── bin1_Lactiplantibacillus_plantarum.consensus.fasta
│   │   ├── bin1_Lactiplantibacillus_plantarum.filtered.vcf
│   │   ├── prokka/
│   │   │   ├── Lactiplantibacillus_plantarum.faa
│   │   │   ├── Lactiplantibacillus_plantarum.gff
│   ├── bin2_Salmonella_enterica/
│   └── ...
├── visualizations/
│   ├── bin1_Lactiplantibacillus_plantarum/
│   │   ├── cog_distribution.png
│   │   ├── top_pathways.png
│   │   ├── pathway_heatmap.png
│   │   ├── pathway_overlap.png
│   │   ├── analysis_results.xlsx
│   ├── bin2_Salmonella_enterica/
│   └── ...
├── clean_R1.fastq
├── clean_R2.fastq
├── fastp.html
├── fastp.json
├── aligned.bam
├── depth.txt
├── mapped.bam
└── pipeline.log
```

## Dependencies
- **Python Packages**: `biopython`
- **Bioinformatics Tools**: `fastp`, `spades.py`, `minimap2`, `samtools`, `metabat2`, `jgi_summarize_bam_contig_depths`, `kraken2`, `abricate`, `prokka`, `diamond`, `snpEff`
- **Databases**: Kraken2 database, eggNOG 5.0 database, organism-specific reference genomes, SnpEff databases

## Troubleshooting
1. **Kraken2 Database Issues**:
   - Ensure `KRAKEN2_DB_PATH` or `--db_path` points to a valid Kraken2 database.
   - Use a comprehensive database (e.g., Standard or PlusPFP) for better taxonomic resolution:
     ```bash
     export KRAKEN2_DB_PATH=/path/to/standard_kraken2_db
     ```

2. **eggNOG Database Issues**:
   - Verify that `EGGNOG_DATA_DIR` is set and contains the eggNOG database:
     ```bash
     ls -l $EGGNOG_DATA_DIR/eggnog.db
     ```
   - If `eggnog-mapper` fails, try downloading the database again or use `--skip-emapper`.

3. **Reference Genome Errors**:
   - Verify that `--ref_dir` contains `reference.fasta` and `reference.gbk`:
     ```bash
     ls -l ref_dir_Lactiplantibacillus_plantarum/
     ```
   - Check `pipeline.log` for bin-specific taxa to select appropriate references.

4. **Empty FASTQ Files**:
   - Check `bins/binID/` for empty FASTQ files:
     ```bash
     find output_dir/bins -name "*_R1.fastq" -empty
     ```
   - If empty, verify Kraken2 outputs and read mapping:
     ```bash
     cat output_dir/kraken2_output/bin1_kraken2_output.txt | grep '^C' | cut -f 3 | sort | uniq
     ```

5. **General Debugging**:
   - Review `pipeline.log` and `metagenomics_pipeline.log` for errors.
   - Use `--verbose` for detailed output and `--dry-run` to preview commands.

## Contributing
Contributions are welcome! Please submit a pull request or open an issue on the repository for bugs, feature requests, or improvements.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.
