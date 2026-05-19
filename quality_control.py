#!/usr/bin/env python3
"""
Quality Control Module for Microbial Genome Assembly Pipeline

This module provides comprehensive quality control capabilities including:
- FastQ quality assessment (fastp)
- Assembly quality evaluation (QUAST)
- Read mapping statistics (samtools)
- MultiQC aggregated reporting

Quality Control Checkpoints:
  QC-1: FastQ Quality Assessment (fastp) - Pre-processing validation
  QC-2: Assembly Quality Evaluation (QUAST) - Consensus sequence assessment
  QC-3: Read Mapping Statistics - Coverage and alignment metrics
  QC-4: Comprehensive Summary Report - JSON + Text formats
  QC-5: MultiQC Aggregated Report - Interactive HTML dashboard
"""

import logging
import os
import subprocess
import json
import sys
from pathlib import Path

def setup_qc_logging(output_dir):
    """Setup logging for QC module."""
    log_file = os.path.join(output_dir, "qc.log")
    handler = logging.FileHandler(log_file)
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(handler)
    return log_file

def run_command_check(cmd, step_name, output_file=None):
    """Run a command and check for success."""
    try:
        if output_file:
            with open(output_file, 'w') as f:
                result = subprocess.run(cmd, shell=True, check=True, text=True, 
                                      stdout=f, stderr=subprocess.PIPE)
        else:
            result = subprocess.run(cmd, shell=True, check=True, text=True, 
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"{step_name}: SUCCESS")
        return True, result.stdout if not output_file else ""
    except subprocess.CalledProcessError as e:
        logging.error(f"{step_name}: FAILED - {e.stderr}")
        return False, e.stderr

def assess_fastq_quality(fastq1, fastq2, output_dir, threads=8, resume=False):
    """
    Assess FastQ quality using fastp.
    
    Returns:
        Tuple: (clean_fastq1_path, clean_fastq2_path, qc_metrics_dict)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    clean_r1 = os.path.join(output_dir, "clean_R1.fastq.gz")
    clean_r2 = os.path.join(output_dir, "clean_R2.fastq.gz")
    fastp_html = os.path.join(output_dir, "fastp.html")
    fastp_json = os.path.join(output_dir, "fastp.json")
    
    if resume and os.path.exists(clean_r1) and os.path.exists(clean_r2):
        logging.info(f"FastQ QC outputs found, skipping fastp run")
        metrics = _parse_fastp_json(fastp_json) if os.path.exists(fastp_json) else {}
        return clean_r1, clean_r2, metrics
    
    logging.info(f"Running fastp on {fastq1} and {fastq2}")
    
    cmd = (
        f"fastp --in1 {fastq1} --in2 {fastq2} "
        f"--out1 {clean_r1} --out2 {clean_r2} "
        f"--html {fastp_html} --json {fastp_json} "
        f"--thread {threads} --adapter_sequence auto "
        f"--cut_front --cut_tail --cut_mean_quality 20 "
        f"--length_required 50 --low_complexity_filter"
    )
    
    success, _ = run_command_check(cmd, "FastQ Quality Assessment (fastp)")
    
    if not success:
        logging.error("fastp failed")
        return None, None, {}
    
    # Parse quality metrics
    metrics = _parse_fastp_json(fastp_json)
    logging.info(f"FastQ QC metrics: {metrics}")
    
    return clean_r1, clean_r2, metrics

def check_assembly_quality_quast(assembly_fasta, reference_fasta, output_dir, threads=8, resume=False):
    """
    Evaluate assembly quality using QUAST.
    
    Returns:
        Dictionary: Assembly quality metrics
    """
    quast_dir = os.path.join(output_dir, "quast_results")
    quast_report = os.path.join(quast_dir, "report.json")
    
    if resume and os.path.exists(quast_report):
        logging.info("QUAST results found, skipping QUAST run")
        metrics = _parse_quast_json(quast_report)
        return metrics
    
    logging.info(f"Running QUAST on {assembly_fasta}")
    
    cmd = (
        f"quast.py {assembly_fasta} -r {reference_fasta} "
        f"-o {quast_dir} -t {threads} --silent --json-output "
        f"--no-icarus --no-snps"
    )
    
    success, _ = run_command_check(cmd, "Assembly Quality Assessment (QUAST)")
    
    if not success:
        logging.error("QUAST failed")
        return {}
    
    metrics = _parse_quast_json(quast_report)
    logging.info(f"Assembly QC metrics: {metrics}")
    
    return metrics

def calculate_bam_statistics(bam_file, output_dir, prefix=""):
    """
    Calculate BAM file statistics using samtools.
    
    Returns:
        Dictionary: Mapping statistics
    """
    os.makedirs(output_dir, exist_ok=True)
    
    stats_file = os.path.join(output_dir, f"{prefix}_samtools.stats")
    flagstat_file = os.path.join(output_dir, f"{prefix}_samtools.flagstat")
    
    logging.info(f"Calculating BAM statistics for {bam_file}")
    
    # Run samtools stats
    cmd_stats = f"samtools stats {bam_file} > {stats_file}"
    success_stats, _ = run_command_check(cmd_stats, "BAM Statistics (samtools stats)")
    
    # Run samtools flagstat
    cmd_flagstat = f"samtools flagstat {bam_file} > {flagstat_file}"
    success_flagstat, _ = run_command_check(cmd_flagstat, "BAM Flagstat (samtools flagstat)")
    
    metrics = {}
    if success_stats:
        metrics = _parse_samtools_stats(stats_file)
    if success_flagstat:
        metrics.update(_parse_samtools_flagstat(flagstat_file))
    
    logging.info(f"Mapping statistics: {metrics}")
    return metrics

def validate_assembly_quality(assembly_metrics, min_n50=5000, max_ngaps=0):
    """
    Validate assembly quality against predefined thresholds.
    
    Args:
        assembly_metrics: Dictionary from QUAST analysis
        min_n50: Minimum N50 value (bp)
        max_ngaps: Maximum number of gaps allowed
    
    Returns:
        Dictionary: Validation results
    """
    validation = {
        "passed": True,
        "warnings": [],
        "errors": []
    }
    
    if not assembly_metrics:
        validation["errors"].append("No assembly metrics available")
        validation["passed"] = False
        return validation
    
    # Check N50
    n50 = assembly_metrics.get("N50", 0)
    if n50 < min_n50:
        validation["warnings"].append(f"Low N50: {n50} bp (threshold: {min_n50} bp)")
    
    # Check number of gaps
    ngaps = assembly_metrics.get("# gaps", 0)
    if ngaps > max_ngaps:
        validation["warnings"].append(f"Gaps detected: {ngaps} (max allowed: {max_ngaps})")
    
    # Check genome size
    genome_size = assembly_metrics.get("Total length", 0)
    if genome_size == 0:
        validation["errors"].append("Genome size is 0")
        validation["passed"] = False
    
    return validation

def generate_qc_summary_report(fastq_metrics, assembly_metrics, mapping_stats, output_file, validation=None):
    """
    Generate comprehensive QC summary report.
    
    Combines FastQ, Assembly, and Mapping QC results into a human-readable report.
    """
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("QUALITY CONTROL SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # FastQ QC
        f.write("1. FASTQ QUALITY ASSESSMENT (fastp)\n")
        f.write("-" * 80 + "\n")
        if fastq_metrics:
            f.write(f"  Read1 Length: {fastq_metrics.get('read1_length', 'N/A')} bp\n")
            f.write(f"  Read2 Length: {fastq_metrics.get('read2_length', 'N/A')} bp\n")
            f.write(f"  Read1 GC: {fastq_metrics.get('read1_gc', 'N/A')}%\n")
            f.write(f"  Read2 GC: {fastq_metrics.get('read2_gc', 'N/A')}%\n")
            f.write(f"  Adapter-trimmed: {fastq_metrics.get('adapter_trimmed', 'N/A')} reads\n")
            f.write(f"  Quality-filtered: {fastq_metrics.get('quality_filtered', 'N/A')} reads\n")
        else:
            f.write("  No metrics available\n")
        f.write("\n")
        
        # Assembly QC
        f.write("2. ASSEMBLY QUALITY ASSESSMENT (QUAST)\n")
        f.write("-" * 80 + "\n")
        if assembly_metrics:
            f.write(f"  Total Length: {assembly_metrics.get('Total length', 'N/A')} bp\n")
            f.write(f"  N50: {assembly_metrics.get('N50', 'N/A')} bp\n")
            f.write(f"  L50: {assembly_metrics.get('L50', 'N/A')}\n")
            f.write(f"  Number of Contigs: {assembly_metrics.get('# contigs', 'N/A')}\n")
            f.write(f"  Number of Gaps: {assembly_metrics.get('# gaps', 'N/A')}\n")
            f.write(f"  GC Content: {assembly_metrics.get('GC (%)', 'N/A')}%\n")
        else:
            f.write("  No metrics available\n")
        f.write("\n")
        
        # Mapping QC
        f.write("3. READ MAPPING STATISTICS (samtools)\n")
        f.write("-" * 80 + "\n")
        if mapping_stats:
            f.write(f"  Total Reads: {mapping_stats.get('total_reads', 'N/A')}\n")
            f.write(f"  Mapped Reads: {mapping_stats.get('mapped_reads', 'N/A')}\n")
            f.write(f"  Mapping Rate: {mapping_stats.get('mapping_rate', 'N/A')}%\n")
            f.write(f"  Mean Coverage: {mapping_stats.get('mean_coverage', 'N/A')}x\n")
            f.write(f"  Duplicates: {mapping_stats.get('duplicates', 'N/A')}%\n")
        else:
            f.write("  No metrics available\n")
        f.write("\n")
        
        # Validation Summary
        f.write("4. QUALITY VALIDATION SUMMARY\n")
        f.write("-" * 80 + "\n")
        if validation:
            status = "✓ PASSED" if validation['passed'] else "✗ FAILED"
            f.write(f"  Status: {status}\n")
            if validation['warnings']:
                f.write("  Warnings:\n")
                for warn in validation['warnings']:
                    f.write(f"    - {warn}\n")
            if validation['errors']:
                f.write("  Errors:\n")
                for err in validation['errors']:
                    f.write(f"    - {err}\n")
        f.write("\n")
        
        f.write("=" * 80 + "\n")
    
    logging.info(f"QC summary report written to {output_file}")

def generate_multiqc_report(results_dir, output_dir, project_name="Pipeline", resume=False):
    """
    Generate MultiQC aggregated report.
    
    Combines all QC tool outputs into a single interactive HTML report.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    multiqc_report = os.path.join(output_dir, "multiqc_report.html")
    
    if resume and os.path.exists(multiqc_report):
        logging.info("MultiQC report found, skipping generation")
        return multiqc_report
    
    logging.info(f"Running MultiQC on {results_dir}")
    
    cmd = (
        f"multiqc {results_dir} -o {output_dir} "
        f"-n multiqc_report --title '{project_name}' "
        f"--force --quiet"
    )
    
    success, _ = run_command_check(cmd, "MultiQC Report Generation")
    
    if not success:
        logging.warning("MultiQC failed, trying alternative approach")
        # Continue even if MultiQC fails
        multiqc_report = None
    
    return multiqc_report

# ===== PARSING FUNCTIONS =====

def _parse_fastp_json(json_file):
    """Parse fastp JSON output."""
    metrics = {}
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            
            summary = data.get('summary', {})
            before = data.get('read1_before_filtering', {})
            after = data.get('read1_after_filtering', {})
            
            metrics['read1_length'] = before.get('mean_length', 0)
            metrics['read1_gc'] = round(before.get('gc_content', 0) * 100, 2)
            metrics['adapter_trimmed'] = summary.get('adapter_trimmed', 0)
            metrics['quality_filtered'] = summary.get('filtering_result', {}).get('low_quality_reads', 0)
            
        logging.info(f"Parsed fastp metrics: {metrics}")
    except Exception as e:
        logging.warning(f"Could not parse fastp JSON: {e}")
    
    return metrics

def _parse_quast_json(json_file):
    """Parse QUAST JSON output."""
    metrics = {}
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            
            # QUAST stores metrics per assembly
            for assembly_data in data.get('results', []):
                for key, value in assembly_data.items():
                    if key not in ['Number of Ns per 100 kbp']:
                        metrics[key] = value
        
        logging.info(f"Parsed QUAST metrics: {metrics}")
    except Exception as e:
        logging.warning(f"Could not parse QUAST JSON: {e}")
    
    return metrics

def _parse_samtools_stats(stats_file):
    """Parse samtools stats output."""
    metrics = {}
    try:
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        key = parts[1].rstrip(':')
                        try:
                            value = float(parts[2])
                            metrics[key] = value
                        except ValueError:
                            metrics[key] = parts[2]
        
        logging.info(f"Parsed samtools stats: {metrics}")
    except Exception as e:
        logging.warning(f"Could not parse samtools stats: {e}")
    
    return metrics

def _parse_samtools_flagstat(flagstat_file):
    """Parse samtools flagstat output."""
    metrics = {}
    try:
        with open(flagstat_file, 'r') as f:
            lines = f.readlines()
            
            if len(lines) > 0:
                parts = lines[0].strip().split()
                metrics['total_reads'] = int(parts[0])
            
            if len(lines) > 3:
                parts = lines[4].strip().split()
                mapped = int(parts[0])
                total = metrics.get('total_reads', 1)
                metrics['mapped_reads'] = mapped
                metrics['mapping_rate'] = round((mapped / total * 100), 2) if total > 0 else 0
            
            if len(lines) > 8:
                parts = lines[8].strip().split()
                duplicates = int(parts[0])
                total = metrics.get('total_reads', 1)
                metrics['duplicates'] = round((duplicates / total * 100), 2) if total > 0 else 0
        
        logging.info(f"Parsed samtools flagstat: {metrics}")
    except Exception as e:
        logging.warning(f"Could not parse samtools flagstat: {e}")
    
    return metrics
