# Quality Control Integration Guide

## Overview

This document describes the comprehensive quality control (QC) system integrated into the Microbial Genome Assembly and Functional Annotation Pipeline. The QC framework ensures that your pipeline produces publication-quality data by validating input reads, assessing assembly quality, and verifying read mapping statistics.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Quality Control Checkpoints](#quality-control-checkpoints)
3. [Quality Thresholds](#quality-thresholds)
4. [Output Files](#output-files)
5. [Report Interpretation](#report-interpretation)
6. [Troubleshooting](#troubleshooting)
7. [Advanced Configuration](#advanced-configuration)

---

## Quick Start

### Enable Quality Control (Default)

```bash
python genomic_pipeline.py \
  -s sample_name \
  -1 input_R1.fastq \
  -2 input_R2.fastq \
  -r reference_directory \
  -d snpeff_db \
  -o output_directory
```

All QC checkpoints run automatically. Results are saved in `output/qc_results/`.

### Skip Quality Control (Testing Only)

```bash
python genomic_pipeline.py \
  -s sample_name \
  -1 input_R1.fastq \
  -2 input_R2.fastq \
  -r reference_directory \
  -d snpeff_db \
  -o output_directory \
  --skip-qc
```

### Resume with Checkpoints

```bash
python genomic_pipeline.py \
  -s sample_name \
  -1 input_R1.fastq \
  -2 input_R2.fastq \
  -r reference_directory \
  -d snpeff_db \
  -o output_directory \
  --resume
```

The pipeline skips completed QC steps and continues from where it left off.

---

## Quality Control Checkpoints

### QC-1: FastQ Quality Assessment (fastp)

**Purpose**: Evaluate and filter input sequencing reads

**Tool**: [fastp](http://opengene.org/fastp/)

**What It Does**:
- Detects and removes adapter sequences (auto-detection)
- Trims low-quality bases (Phred score < 20)
- Removes reads shorter than 50 bp
- Filters low-complexity sequences
- Generates quality statistics

**Output Files**:
```
qc_results/
├── clean_R1.fastq.gz          # Filtered forward reads
├── clean_R2.fastq.gz          # Filtered reverse reads
├── fastp.html                 # Interactive quality report
└── fastp.json                 # Detailed metrics (JSON)
```

**Key Metrics**:
- **Read Length**: Mean length of sequencing reads
- **GC Content**: Percentage of G+C nucleotides
- **Adapter Trimmed**: Number of reads with adapters removed
- **Quality Filtered**: Number of reads removed due to low quality

**Interpretation**:
```
✓ Good quality: 
  - GC content: 40-60% (species-dependent)
  - Adapter trimming: < 5% of reads
  - Quality filtering: < 10% of reads
  
⚠ Warning:
  - GC content: 30-40% or 60-70% (possible contamination)
  - Adapter trimming: 5-15% of reads
  
✗ Poor quality:
  - GC content: < 30% or > 70%
  - Adapter trimming: > 15% of reads
  - Quality filtering: > 25% of reads
```

---

### QC-2: Assembly Quality Assessment (QUAST)

**Purpose**: Evaluate the quality of the consensus/assembled sequence

**Tool**: [QUAST](http://quast.sourceforge.net/)

**What It Does**:
- Calculates N50 and L50 values
- Counts contigs and gaps
- Determines GC content
- Compares against reference genome
- Identifies misassemblies

**Output Files**:
```
qc_results/quast_results/
├── report.html                # Interactive QUAST report
├── report.json                # Metrics in JSON format
├── report.txt                 # Text summary
└── contigs_report/            # Per-contig statistics
```

**Key Metrics**:

| Metric | Description | Good Range | Warning | Poor |
|--------|-------------|------------|---------|------|
| **N50** | Contig length median | > 10 kb | 5-10 kb | < 5 kb |
| **L50** | Number of contigs >= N50 | < 10 | 10-20 | > 20 |
| **Total Length** | Genome assembly size | Near reference | ±5% ref | ±10% ref |
| **# Contigs** | Number of assembled contigs | < 50 | 50-100 | > 100 |
| **# Gaps** | Number of assembly gaps | 0 | 1-5 | > 5 |
| **GC (%)** | GC content | Species-dependent | ±2% species avg | ±5% |

**Interpretation**:

```
Excellent Assembly:
  - N50 > 100 kb
  - L50 < 5
  - Total length near reference (±2%)
  - Contigs < 10
  - No gaps
  - GC content matches reference

Good Assembly:
  - N50 > 50 kb
  - L50 < 10
  - Total length near reference (±5%)
  - Contigs < 50
  - Gaps < 3
  - GC content within ±2% of reference

Acceptable Assembly:
  - N50 > 10 kb
  - L50 < 20
  - Total length near reference (±10%)
  - Contigs < 100
  - Gaps < 10
  - GC content within ±5% of reference

Poor Assembly:
  - N50 < 10 kb
  - L50 > 20
  - Total length differs > ±10% from reference
  - Contigs > 100
  - Gaps > 10
  - GC content differs > ±5% from reference
```

---

### QC-3: Read Mapping Statistics

**Purpose**: Verify that sequencing reads map correctly to the assembly

**Tool**: [samtools](http://www.htslib.org/)

**What It Does**:
- Calculates total number of aligned reads
- Determines mapping rate (% reads successfully mapped)
- Computes mean coverage depth
- Identifies duplicate reads
- Provides alignment quality statistics

**Output Files**:
```
qc_results/
├── genomic_samtools.stats    # Detailed alignment statistics
└── genomic_samtools.flagstat # Read mapping summary
```

**Key Metrics**:

| Metric | Description | Good Range | Warning | Poor |
|--------|-------------|------------|---------|------|
| **Total Reads** | Number of sequencing reads | > 1M | 100K-1M | < 100K |
| **Mapped Reads** | Reads aligning to assembly | > 95% | 85-95% | < 85% |
| **Mean Coverage** | Average sequencing depth | > 30x | 15-30x | < 15x |
| **Duplicates** | PCR duplicate percentage | < 10% | 10-20% | > 20% |

**Interpretation**:

```
Excellent Coverage:
  - Mapping rate > 98%
  - Mean coverage > 100x
  - Duplicates < 5%
  → High-confidence assembly, suitable for variant calling

Good Coverage:
  - Mapping rate > 95%
  - Mean coverage > 30x
  - Duplicates < 10%
  → Sufficient for most analyses

Acceptable Coverage:
  - Mapping rate > 85%
  - Mean coverage > 15x
  - Duplicates < 20%
  → May be suitable, use with caution

Poor Coverage:
  - Mapping rate < 85%
  - Mean coverage < 15x
  - Duplicates > 20%
  → Not recommended for variant calling
  → Consider re-sequencing or quality control
```

---

### QC-4: Comprehensive Summary Report

**Purpose**: Provide human-readable summary of all QC metrics

**Format**: Plain text

**Output Files**:
```
qc_results/qc_summary_report.txt
```

**Content**:
- FastQ QC summary
- Assembly QC metrics
- Mapping statistics
- Quality validation status (PASS/FAIL)
- Warnings and errors

**Example Report**:
```
================================================================================
QUALITY CONTROL SUMMARY REPORT
================================================================================

1. FASTQ QUALITY ASSESSMENT (fastp)
--------------------------------------------------------------------------------
  Read1 Length: 150.0 bp
  Read1 GC: 52.3%
  Adapter-trimmed: 125,000 reads
  Quality-filtered: 85,000 reads

2. ASSEMBLY QUALITY ASSESSMENT (QUAST)
--------------------------------------------------------------------------------
  Total Length: 4,850,000 bp
  N50: 125,000 bp
  L50: 15
  Number of Contigs: 42
  Number of Gaps: 2
  GC Content: 50.5%

3. READ MAPPING STATISTICS (samtools)
--------------------------------------------------------------------------------
  Total Reads: 5,000,000
  Mapped Reads: 4,900,000
  Mapping Rate: 98.0%
  Mean Coverage: 52.3x
  Duplicates: 3.5%

4. QUALITY VALIDATION SUMMARY
--------------------------------------------------------------------------------
  Status: ✓ PASSED
  Warnings: None
  Errors: None

================================================================================
```

---

### QC-5: MultiQC Aggregated Report

**Purpose**: Generate interactive HTML dashboard combining all QC tools

**Tool**: [MultiQC](https://multiqc.info/)

**Features**:
- Integrates fastp, QUAST, samtools results
- Interactive visualizations
- Comparative analysis across samples
- Professional publication-ready formatting

**Output Files**:
```
qc_results/
├── multiqc_report.html        # Main interactive report
├── multiqc_data/              # Report data (JSON, TSV)
└── multiqc.log                # MultiQC execution log
```

**How to View**:
1. Open `qc_results/multiqc_report.html` in a web browser
2. Navigate between sections using the sidebar
3. Hover over plots for detailed information
4. Download plots as images using the camera icon

**Sections**:
- **General Statistics**: Summary table of key metrics
- **FastP**: Read quality visualizations
- **QUAST**: Assembly quality plots
- **Samtools**: Coverage and mapping statistics
- **Custom Metrics**: Pipeline-specific information

---

## Quality Thresholds

The pipeline uses the following quality thresholds to determine if data passes QC:

### Assembly Quality Thresholds

```python
N50 minimum:           5,000 bp
Maximum gaps:          0 (no gaps allowed)
GC content tolerance:  ±5% from species average
```

### Mapping Quality Thresholds

```python
Minimum mapping rate:  95%
Minimum mean coverage: 30x
Maximum duplicates:    20%
```

### FastQ Quality Thresholds

```python
Minimum read length:   50 bp (after trimming)
Quality threshold:     Phred 20 (Q20)
GC content tolerance:  ±15% from expected
```

### Custom Thresholds

To modify thresholds, edit `quality_control.py`:

```python
# In validate_assembly_quality() function:
def validate_assembly_quality(assembly_metrics, min_n50=5000, max_ngaps=0):
    """Customize thresholds here"""
    pass
```

Then modify your call:
```python
validation = validate_assembly_quality(
    assembly_metrics, 
    min_n50=10000,      # Custom N50 threshold
    max_ngaps=2         # Allow up to 2 gaps
)
```

---

## Output Files

### QC Results Directory Structure

```
output/
├── qc_results/                          # All QC outputs
│   ├── clean_R1.fastq.gz                # QC-1: Filtered reads
│   ├── clean_R2.fastq.gz
│   ├── fastp.html                       # QC-1: Interactive report
│   ├── fastp.json                       # QC-1: Metrics (JSON)
│   ├── quast_results/                   # QC-2: Assembly assessment
│   │   ├── report.html
│   │   ├── report.json
│   │   ├── report.txt
│   │   └── contigs_report/
│   ├── genomic_samtools.stats           # QC-3: Detailed statistics
│   ├── genomic_samtools.flagstat        # QC-3: Mapping summary
│   ├── qc_summary_report.txt            # QC-4: Text summary
│   ├── multiqc_report.html              # QC-5: Interactive dashboard
│   ├── multiqc_data/                    # QC-5: Report data
│   └── qc.log                           # QC execution log
│
├── sample_R1.paired.fastq               # Paired reads
├── sample_R2.paired.fastq
├── sample_R1.sorted.fastq               # Sorted paired reads
├── sample_R2.sorted.fastq
├── sample.sorted.bam                    # Alignment file
├── sample.consensus.fasta               # Final consensus sequence
├── sample.filtered.vcf.gz                # Filtered variants
├── sample.annotated.vcf                 # Annotated variants
├── sample_snpEff_summary.html           # SnpEff report
├── sample_snpEff_summary.csv            # SnpEff metrics
├── prokka/                              # Gene annotation results
│   ├── sample.gff
│   ├── sample.gbk
│   ├── sample.faa
│   └── ...
└── [other outputs]
```

### File Descriptions

| File | Source | Content | Size |
|------|--------|---------|------|
| `clean_R1.fastq.gz` | fastp | Quality-filtered reads (gzipped) | Variable |
| `fastp.html` | fastp | Interactive quality report | 1-5 MB |
| `report.html` (QUAST) | QUAST | Assembly assessment with plots | 2-10 MB |
| `genomic_samtools.stats` | samtools | Alignment statistics | 10-100 KB |
| `qc_summary_report.txt` | Custom | Human-readable summary | 5-20 KB |
| `multiqc_report.html` | MultiQC | Aggregated dashboard | 5-20 MB |

---

## Report Interpretation

### Interpreting fastp Report

1. **Open `fastp.html`** in a web browser
2. **Quality Curves**: 
   - X-axis: Position in read
   - Y-axis: Quality score
   - Green: Good quality (Q ≥ 30)
   - Yellow: Moderate (Q 20-30)
   - Red: Poor (Q < 20)
3. **Adapter Content**: Should be < 5% for good data
4. **Duplication Rate**: Should be < 30% (depends on organism)

**What to Look For**:
```
✓ Good Signs:
  - Quality lines remain high across entire read
  - Minimal adapter contamination
  - Even duplication rates

✗ Bad Signs:
  - Quality drops sharply at read ends
  - High adapter content (> 20%)
  - Extremely high duplication (> 50%)
```

### Interpreting QUAST Report

1. **Open `report.html`** in a web browser
2. **Summary Table**: Shows key metrics
3. **Contigs**: Sorted by length
4. **GC Content**: Compared to reference
5. **Icarus**: Visual contig alignment (if reference provided)

**What to Look For**:
```
✓ Good Signs:
  - N50 > 50 kb
  - Few contigs (< 50)
  - GC content matches reference
  - Minimal gaps
  - High NG50 (if reference available)

✗ Bad Signs:
  - N50 < 10 kb (fragmented)
  - Many contigs (> 100)
  - GC content differs significantly from reference
  - Multiple gaps (> 10)
  - Low NG50 relative to N50
```

### Interpreting MultiQC Report

1. **General Statistics Tab**: Overview of all samples
2. **FastP Section**: Quality before/after filtering
3. **QUAST Section**: Assembly metrics and plots
4. **Samtools Section**: Mapping coverage profiles
5. **Custom Content**: Pipeline-specific metrics

**Navigation Tips**:
- Use sidebar to jump between sections
- Hover over plots for details
- Click on column headers to sort tables
- Download data using the download button

---

## Troubleshooting

### Problem: FastQ QC Shows Poor Quality

**Symptoms**:
```
- GC content < 30% or > 70%
- > 25% quality filtering
- Uneven quality across reads
```

**Solutions**:
1. Check sequencing machine logs for errors
2. Verify sample preparation protocol
3. Consider re-sequencing with QC
4. Check for contamination:
   ```bash
   kraken2 --db kraken_db clean_R1.fastq.gz | head -100
   ```

---

### Problem: Low Assembly Quality (N50 < 10 kb)

**Symptoms**:
```
- N50 < 10 kb
- High number of contigs (> 100)
- Multiple gaps
```

**Solutions**:
1. Check input read depth:
   ```bash
   zcat clean_R1.fastq.gz | wc -l
   # Divide by 4 to get number of reads
   ```
2. Increase sequencing depth (usually needs > 30x coverage)
3. Verify reference genome is correct
4. Try different assembly parameters in SPAdes/Bowtie2

---

### Problem: Low Mapping Rate (< 85%)

**Symptoms**:
```
- Mapping rate < 85%
- Mean coverage < 15x
- Many unmapped reads
```

**Solutions**:
1. Verify reference genome is correct:
   ```bash
   samtools flagstat sample.sorted.bam
   ```
2. Check for contamination in input data
3. Verify read format (FASTQ vs FASTA)
4. Check if reference needs indexing:
   ```bash
   samtools faidx reference.fasta
   bowtie2-build reference.fasta reference_index
   ```

---

### Problem: MultiQC Report Not Generated

**Symptoms**:
```
- multiqc_report.html not created
- Error in pipeline log
```

**Solutions**:
1. Ensure MultiQC is installed:
   ```bash
   conda install -c bioconda multiqc
   multiqc --version
   ```
2. Check directory permissions:
   ```bash
   ls -la output/qc_results/
   ```
3. Run MultiQC manually:
   ```bash
   multiqc output/qc_results/ -o output/qc_results/
   ```
4. Check for MultiQC errors:
   ```bash
   multiqc output/qc_results/ -v  # Verbose mode
   ```

---

### Problem: Disk Space Issues During QC

**Symptoms**:
```
- Pipeline stops with "No space left on device"
- Clean reads files are very large
```

**Solutions**:
1. Clean temporary files:
   ```bash
   rm -f output/*sorted.fastq
   rm -f output/*.sam
   ```
2. Compress clean reads:
   ```bash
   gzip output/qc_results/clean_*.fastq
   ```
3. Enable temp file cleanup:
   ```bash
   python genomic_pipeline.py ... # (enabled by default)
   ```
4. Check disk space:
   ```bash
   df -h
   ```

---

## Advanced Configuration

### Custom Quality Control Functions

You can modify the QC module for custom thresholds:

```python
# In genomic_pipeline.py:
from quality_control import validate_assembly_quality

# Custom validation
assembly_validation = validate_assembly_quality(
    assembly_metrics,
    min_n50=10000,      # Increase N50 threshold
    max_ngaps=2         # Allow up to 2 gaps
)
```

### Skip Individual QC Checkpoints

Modify the `main()` function in `genomic_pipeline.py`:

```python
# QC-1: Skip FastQ assessment
# Comment out this block:
# logging.info("[QC-1] Assessing FastQ quality...")
# clean_fastq1, clean_fastq2, fastq_qc = assess_fastq_quality(...)

# QC-2: Skip Assembly assessment
# Comment out:
# assembly_metrics = check_assembly_quality_quast(...)

# etc.
```

### Custom Report Template

Modify `generate_qc_summary_report()` in `quality_control.py`:

```python
def generate_qc_summary_report(fastq_metrics, assembly_metrics, 
                              mapping_stats, output_file, validation=None):
    with open(output_file, 'w') as f:
        # Add custom report sections here
        f.write("CUSTOM QC REPORT\n")
        # Add your custom metrics
```

### Integrate Additional QC Tools

Add new tools to the QC module:

```python
def assess_contamination(fastq_file, output_dir, db_path):
    """Assess contamination using Kraken2"""
    kraken_output = os.path.join(output_dir, "kraken2.txt")
    cmd = f"kraken2 --db {db_path} {fastq_file} > {kraken_output}"
    run_command_check(cmd, "Contamination Assessment (Kraken2)")
    return kraken_output

# Call from genomic_pipeline.py:
contamination_report = assess_contamination(clean_fastq1, qc_dir, kraken_db)
```

---

## References

### Tool Documentation

- **fastp**: http://opengene.org/fastp/ - Quality control and preprocessing
- **QUAST**: http://quast.sourceforge.net/ - Assembly quality assessment
- **samtools**: http://www.htslib.org/ - BAM/SAM file processing
- **MultiQC**: https://multiqc.info/ - Report aggregation

### Recommended Reading

- [Illumina Quality Scores](https://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/understanding-illumina-quality-scores-technical-note-370-2017-004.pdf)
- [Genome Assembly Best Practices](https://github.com/PacificBiosciences/DevNet/wiki/Genome-Assembly-Best-Practices)
- [Variant Calling Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Variant-Calling-Best-Practices)

### Citation

If you use this QC framework in your research, please cite:

```
Tools used in this pipeline:
- fastp: Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast 
  all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.
  
- QUAST: Gurevich, A., Saveliev, V., Vyatkin, A., & Testov, V. (2013). 
  QUAST: quality assessment tool for genome assemblies. Bioinformatics, 
  29(8), 1072-1075.
  
- MultiQC: Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). 
  MultiQC: summarize analysis results for multiple tools and samples in a 
  single report. Bioinformatics, 32(19), 3047-3048.
```

---

## Support

For issues, questions, or suggestions:

1. Check this documentation
2. Review pipeline logs: `output/qc_results/qc.log`
3. Run in verbose mode: `--verbose` flag
4. Check tool documentation (links above)
5. Open an issue on GitHub with logs and detailed description

---

**Last Updated**: 2026-05-19
**Version**: 1.0
**Maintainer**: Rajindra04
