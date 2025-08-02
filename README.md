# Metagenomics Analysis Pipeline

This repository contains a comprehensive pipeline for metagenomics analysis, variant calling, and visualization of genomic annotations. The pipeline integrates three main components:

1. **Metagenomics Pipeline**: Processes paired-end FASTQ files to perform quality control, de novo assembly, binning, taxonomic classification, and antimicrobial resistance (AMR) detection.
2. **Genomic Pipeline**: Performs variant calling and functional annotation on segregated FASTQ files from the metagenomics pipeline.
3. **Visualization Pipeline**: Generates visualizations of functional annotations, including COG category distributions, metabolic pathway analyses, and gene-pathway associations.

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Running the Full Pipeline](#running-the-full-pipeline)
  - [Running Individual Pipelines](#running-individual-pipelines)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [Contributing](#contributing)
- [License](#license)

## Features

- **Quality Control**: Uses `fastp` to clean paired-end FASTQ files.
- **De Novo Assembly**: Assembles contigs using `SPAdes` in metagenomic mode.
- **Binning**: Groups contigs into bins using `MetaBAT2`.
- **Taxonomic Classification**: Assigns taxonomy to bins using `Kraken2`.
- **AMR and Virulence Detection**: Identifies AMR, virulence, and plasmid genes with `ABRicate`.
- **Variant Calling**: Maps reads to a reference genome, calls variants with `bcftools`, and annotates with `SnpEff`.
- **Functional Annotation**: Annotates genes with `Prokka` and `eggNOG-mapper`.
- **Visualization**: Generates plots for COG distributions, pathway analyses, and gene-pathway associations using `matplotlib` and `seaborn`.
- **Resume Capability**: Supports resuming from completed steps to avoid redundant processing.
- **Modular Design**: Allows running individual pipelines (metagenomics, genomic, visualization) or the full workflow.

## Requirements

- **Operating System**: Linux or macOS (Windows not supported due to dependency compatibility).
- **Conda**: Required for environment management.
- **Kraken2 Database**: Must be specified via `--db_path` or the `KRAKEN2_DB_PATH` environment variable.
- **Reference Files**: Required for the genomic pipeline (`reference.fasta` and `reference.gbk` in `--ref_dir`).

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

4. **Set Up Kraken2 Database** (if not already available):
   - Specify the path to an existing Kraken2 database via the `KRAKEN2_DB_PATH` environment variable:
     ```bash
     export KRAKEN2_DB_PATH=/path/to/kraken2/database
     ```
   - Alternatively, the pipeline can download the MiniKraken2 database automatically if not found (requires internet access).

## Usage

The pipeline is controlled via `main.py`, which allows running the full pipeline or individual components.

### Running the Full Pipeline

To run the entire pipeline (metagenomics, genomic, and visualization):
```bash
python main.py --pipeline all \
    --fastq1 input_R1.fastq \
    --fastq2 input_R2.fastq \
    --output output_dir \
    --ref_dir ref_dir \
    --snpeff_db mydb \
    --threads 8 \
    --resume
```

**Arguments**:
- `--pipeline all`: Runs all pipelines sequentially.
- `--fastq1` and `--fastq2`: Paths to paired-end FASTQ files.
- `--output`: Output directory for results.
- `--ref_dir`: Directory containing `reference.fasta` and `reference.gbk`.
- `--snpeff_db`: SnpEff database name.
- `--threads`: Number of CPU threads (default: 8).
- `--resume`: Resume from completed steps.
- `--db_path`: Path to Kraken2 database (optional, overrides `KRAKEN2_DB_PATH`).
- `--organism`: Organism code for KEGG (default: `eco`).
- `--keep-temp`: Keep temporary files.
- `--dry-run`: Print commands without executing.
- `--verbose`: Enable verbose logging.
- `--skip-emapper`: Skip eggNOG-mapper if the database is unavailable.

### Running Individual Pipelines

- **Metagenomics Pipeline**:
  ```bash
  python main.py --pipeline metagenomics \
      --fastq1 input_R1.fastq \
      --fastq2 input_R2.fastq \
      --output output_dir \
      --threads 8 \
      --resume
  ```

- **Genomic Pipeline** (uses segregated FASTQ files from metagenomics output):
  ```bash
  python main.py --pipeline genomic \
      --output output_dir \
      --ref_dir ref_dir \
      --snpeff_db mydb \
      --threads 8 \
      --resume
  ```

- **Visualization Pipeline** (uses annotation files from genomic output):
  ```bash
  python main.py --pipeline visualize \
      --output output_dir
  ```

## Directory Structure

The pipeline generates the following output structure in the specified `--output` directory:
```
output_dir/
├── bins/                       # Metagenomic bins (bin*.fa)
├── kraken2_output/             # Kraken2 taxonomic classification results
├── spades_output/              # SPAdes assembly results (contigs.fasta)
├── *.fastq                    # Segregated FASTQ files (e.g., Escherichia_coli_R*.fastq)
├── genomic/
│   ├── organism1/             # Genomic pipeline results for organism1
│   │   ├── organism1.annotations.txt
│   │   ├── organism1.consensus.fasta
│   │   ├── organism1.filtered.vcf
│   │   ├── prokka/
│   │   └── ...
│   ├── organism2/
│   └── ...
├── visualizations/
│   ├── organism1/             # Visualization results for organism1
│   │   ├── cog_distribution.png
│   │   ├── top_pathways.png
│   │   ├── pathway_heatmap.png
│   │   ├── pathway_overlap.png
│   │   └── analysis_results.xlsx
│   ├── organism2/
│   └── ...
├── fastp.html                 # Fastp quality control report
├── fastp.json
└── pipeline.log               # Pipeline execution log
```

## Dependencies

The pipeline relies on the following tools and libraries, managed via `environment.yml`:

- **Bioconda Tools**:
  - fastp (0.23.4)
  - spades (4.0.0)
  - minimap2 (2.28)
  - samtools (1.21)
  - metabat2 (2.15)
  - kraken2 (2.1.3)
  - abricate (1.0.0)
  - bowtie2 (2.5.4)
  - bcftools (1.16)
  - snpeff (5.1)
  - prokka (1.14.6)
  - eggnog-mapper (2.1.9)

- **Python Libraries** (via pip):
  - argparse
  - staramr
  - bioservices
  - matplotlib
  - pandas
  - seaborn
  - xlsxwriter
  - matplotlib-venn

- **Other**:
  - r-base (4.2.3)
  - rpy2 (3.5.11)

See `environment.yml` for the complete list.

## Contributing

Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature/your-feature`).
3. Commit your changes (`git commit -m 'Add your feature'`).
4. Push to the branch (`git push origin feature/your-feature`).
5. Open a pull request.

Please ensure your code follows the existing style and includes appropriate tests.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
