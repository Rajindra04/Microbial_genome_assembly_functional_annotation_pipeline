import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from bioservices import KEGG
import time
import os
from matplotlib_venn import venn2, venn3
from itertools import combinations
from time import sleep
from random import uniform
import argparse

# Set style for plots
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 150

def parse_annotations(file_path):
    """Parse the annotations file and clean KEGG column"""
    print(f"Reading file: {file_path}")
    data = []
    invalid_kos = set()
    
    try:
        with open(file_path, 'r') as f:
            header = None
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    if line.startswith('#') and not header:
                        header = line  # Capture the header for reference
                    continue
                parts = line.split('\t')
                # Ensure the line has enough columns (at least 7 for KEGG)
                if len(parts) < 7:
                    print(f"Warning: Skipping line with insufficient columns: {line}")
                    continue
                entry = {
                    'Gene': parts[0],
                    'COG': parts[1] if len(parts) > 1 else '-',
                    'Description': parts[2] if len(parts) > 2 else '',
                    'Preferred_name': parts[3] if len(parts) > 3 else '',
                    'GO': parts[4] if len(parts) > 4 else '',
                    'EC': parts[5] if len(parts) > 5 else '',
                    'KEGG': parts[6] if len(parts) > 6 else '',
                    'KEGG_Pathway': parts[7] if len(parts) > 7 else '',
                    'KEGG_Module': parts[8] if len(parts) > 8 else '',
                    'KEGG_Reaction': parts[9] if len(parts) > 9 else '',
                    'KEGG_rclass': parts[10] if len(parts) > 10 else ''
                }
                # Clean KEGG column: move pathway IDs to KEGG_Pathway
                if entry['KEGG'] and entry['KEGG'] != '-':
                    ko_ids = [ko.strip() for ko in entry['KEGG'].split(',') if ko.strip()]
                    valid_kos = []
                    pathway_ids = []
                    for ko in ko_ids:
                        if ko.startswith('ko:K') and ko[4:].isdigit():
                            valid_kos.append(ko)
                        elif (ko.startswith('ko') or ko.startswith('map')) and ko[3:].isdigit():
                            pathway_ids.append(ko)
                            invalid_kos.add(ko)
                        else:
                            invalid_kos.add(ko)
                    entry['KEGG'] = ','.join(valid_kos) if valid_kos else '-'
                    if pathway_ids:
                        existing_pathways = entry['KEGG_Pathway'].split(',') if entry['KEGG_Pathway'] and entry['KEGG_Pathway'] != '-' else []
                        entry['KEGG_Pathway'] = ','.join(list(set(existing_pathways + pathway_ids)))
                data.append(entry)
    except FileNotFoundError:
        print(f"Error: Annotations file {file_path} not found")
        raise
    except Exception as e:
        print(f"Error reading annotations file {file_path}: {e}")
        raise
    
    df = pd.DataFrame(data)
    print(f"Parsed {len(df)} genes")
    if invalid_kos:
        print(f"Warning: Found invalid KO numbers in KEGG column, moved to KEGG_Pathway or skipped: {', '.join(sorted(invalid_kos))}")
    if len(df) == 0:
        print(f"Warning: No valid gene entries found in {file_path}")
    return df

def plot_cog_distribution(df, output_dir):
    """Plot and save COG category distribution with proper legend"""
    if df.empty:
        print("Warning: Empty DataFrame, skipping COG distribution plot")
        return pd.DataFrame()  # Return empty DataFrame to avoid further errors
    
    print("Counting COG categories...")
    cog_counts = defaultdict(int)
    for cogs in df['COG']:
        for cog in str(cogs).split(','):  # Ensure COG is split correctly
            cog = cog.strip()
            if cog:  # Only count non-empty COGs
                cog_counts[cog] += 1
    
    cog_descriptions = {
        'A': 'RNA processing and modification',
        'B': 'Chromatin structure and dynamics',
        'C': 'Energy production and conversion',
        'D': 'Cell cycle control, cell division',
        'E': 'Amino acid transport and metabolism',
        'F': 'Nucleotide transport and metabolism',
        'G': 'Carbohydrate transport and metabolism',
        'H': 'Coenzyme transport and metabolism',
        'I': 'Lipid transport and metabolism',
        'J': 'Translation, ribosomal structure',
        'K': 'Transcription',
        'L': 'Replication, recombination and repair',
        'M': 'Cell wall/membrane/envelope biogenesis',
        'N': 'Cell motility',
        'O': 'Posttranslational modification',
        'P': 'Inorganic ion transport and metabolism',
        'Q': 'Secondary metabolites biosynthesis',
        'R': 'General function prediction',
        'S': 'Function unknown',
        'T': 'Signal transduction mechanisms',
        'U': 'Intracellular trafficking',
        'V': 'Defense mechanisms',
        '-': 'Not assigned'
    }
    
    cog_data = []
    for cog, count in cog_counts.items():
        cog_data.append({
            'COG': cog,
            'Count': count,
            'Percentage': 100 * count / len(df),
            'Description': cog_descriptions.get(cog, 'Unknown')
        })
    
    cog_df = pd.DataFrame(cog_data).sort_values('Count', ascending=False)
    print(f"Found {len(cog_df)} COG categories")
    
    plt.figure(figsize=(20, 10))  # Increased width for 54 categories
    # Plot top 20 categories to avoid overcrowding, keep full legend
    ax = sns.barplot(x='COG', y='Count', hue='COG', data=cog_df.head(20), palette='viridis', dodge=True)
    
    for p in ax.patches:
        height = p.get_height()
        ax.annotate(f"{height}\n({height/len(df)*100:.1f}%)", 
                   (p.get_x() + p.get_width()/2., height),
                   ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')
    
    # Create custom legend with colors from the viridis palette for all 54 categories
    colors = sns.color_palette('viridis', len(cog_df))
    handles = [plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=colors[i], 
                          markersize=10, label=f"{cog_df['COG'].iloc[i]}: {cog_df['Description'].iloc[i]}") 
               for i in range(len(cog_df))]
    ax.legend(handles=handles, title='COG Categories', bbox_to_anchor=(1.05, 1), loc='upper left', 
              frameon=True, fontsize='small')
    
    plt.title('COG Functional Categories Distribution (Top 20)', pad=20)
    plt.xlabel('COG Category')
    plt.ylabel('Number of Genes')
    plt.xticks(rotation=90, fontsize=8)  # Rotate and reduce font size for readability
    plt.tight_layout(pad=2.0)
    cog_plot_path = os.path.join(output_dir, 'cog_distribution.png')
    plt.savefig(cog_plot_path, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"COG distribution plot saved to {cog_plot_path}")
    
    return cog_df

def get_kegg_pathway_names(pathway_ids, output_dir, cache_file='kegg_cache.csv', max_retries=3):
    """Get pathway names with caching and retry logic"""
    cache_file = os.path.join(output_dir, cache_file)  # Ensure cache file is in output_dir
    cache = {}
    if os.path.exists(cache_file):
        try:
            cache_df = pd.read_csv(cache_file)
            cache = dict(zip(cache_df['pathway_id'], cache_df['pathway_name']))
            print(f"Loaded {len(cache)} pathway names from cache")
        except Exception as e:
            print(f"Error reading cache file: {e}")
    
    k = KEGG()
    pathway_names = {}
    new_entries = {}
    
    for idx, pathway_id in enumerate(pathway_ids):
        if idx % 10 == 0:
            print(f"Fetching pathway names: {idx + 1}/{len(pathway_ids)}")
        if pathway_id in cache:
            pathway_names[pathway_id] = cache[pathway_id]
            continue
        
        if not (pathway_id.startswith('ko') or pathway_id.startswith('map')) or not pathway_id[3:].isdigit():
            pathway_names[pathway_id] = f"Invalid ID ({pathway_id})"
            new_entries[pathway_id] = f"Invalid ID ({pathway_id})"
            continue
        
        print(f"Fetching name for pathway {pathway_id}")
        for attempt in range(max_retries):
            try:
                pathway_info = k.get(pathway_id)
                if pathway_info and isinstance(pathway_info, str):
                    name_lines = [line for line in pathway_info.split('\n') if line.startswith('NAME')]
                    if name_lines:
                        name = name_lines[0].replace('NAME', '').strip()
                        pathway_names[pathway_id] = name
                        new_entries[pathway_id] = name
                    else:
                        pathway_names[pathway_id] = f"No name found ({pathway_id})"
                        new_entries[pathway_id] = f"No name found ({pathway_id})"
                else:
                    pathway_names[pathway_id] = f"No data returned ({pathway_id})"
                    new_entries[pathway_id] = f"No data returned ({pathway_id})"
                sleep(1)
                break
            except Exception as e:
                print(f"Attempt {attempt + 1} failed for {pathway_id}: {str(e)}")
                if attempt == max_retries - 1:
                    pathway_names[pathway_id] = f"Error ({pathway_id}: {str(e)})"
                    new_entries[pathway_id] = f"Error ({pathway_id}: {str(e)})"
                sleep(2 ** attempt + uniform(0, 0.1))
    
    if new_entries:
        try:
            new_cache = pd.DataFrame({
                'pathway_id': list(new_entries.keys()),
                'pathway_name': list(new_entries.values())
            })
            if os.path.exists(cache_file):
                existing_cache = pd.read_csv(cache_file)
                new_cache = pd.concat([existing_cache, new_cache]).drop_duplicates('pathway_id')
            new_cache.to_csv(cache_file, index=False)
            print(f"Cached {len(new_entries)} new pathway names to {cache_file}")
        except Exception as e:
            print(f"Error writing to cache file: {e}")
    
    return pathway_names

def analyze_pathways(df, output_dir, use_ko_mapping=True):
    """Comprehensive pathway analysis with visualizations, mimicking KEGG Mapper Reconstruct"""
    if df.empty:
        print("Warning: Empty DataFrame, skipping pathway analysis")
        return pd.DataFrame(), {}, pd.DataFrame()  # Return empty results to avoid further errors
    
    print("Starting pathway analysis...")
    pathway_genes = defaultdict(list)
    gene_pathways = defaultdict(list)
    k = KEGG()
    
    # Cache for KO-to-pathway mappings
    ko_cache_file = os.path.join(output_dir, 'ko_pathway_cache.csv')  # Ensure cache file is in output_dir
    ko_cache = {}
    if os.path.exists(ko_cache_file):
        try:
            ko_cache_df = pd.read_csv(ko_cache_file)
            ko_cache = dict(zip(ko_cache_df['ko_id'], ko_cache_df['pathways'].apply(lambda x: x.split(',') if x else [])))
            print(f"Loaded {len(ko_cache)} KO-to-pathway mappings from cache")
        except Exception as e:
            print(f"Error reading KO cache file: {e}")
    
    new_ko_entries = {}
    
    # Convert KO numbers to KEGG Mapper format (KXXXXX)
    ko_map = {}
    for ko in df['KEGG'].str.split(',').explode().str.strip().unique():
        if ko and ko != '-' and ko.startswith('ko:K') and ko[4:].isdigit():
            ko_map[ko] = 'K' + ko[4:]
    
    print(f"Found {len(ko_map)} unique valid KO numbers")
    
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"Processing gene {idx + 1}/{len(df)}")
        pathways = set()
        # Process KEGG_Pathway column
        if row['KEGG_Pathway'] and row['KEGG_Pathway'] != '-':
            pathways.update([p.strip() for p in row['KEGG_Pathway'].split(',') if p.strip()])
        
        # Process KEGG (KO) column if enabled
        if use_ko_mapping and row['KEGG'] and row['KEGG'] != '-':
            ko_ids = [ko.strip() for ko in row['KEGG'].split(',') if ko.strip()]
            for ko_id in ko_ids:
                if ko_id not in ko_map:
                    print(f"Skipping invalid KO {ko_id}")
                    continue
                k_number = ko_map[ko_id]
                if ko_id in ko_cache:
                    pathways.update(ko_cache[ko_id])
                    continue
                print(f"Fetching pathways for KO {k_number}")
                try:
                    linked_pathways = k.link('pathway', k_number)
                    ko_pathways = []
                    if linked_pathways:
                        for line in linked_pathways.split('\n'):
                            if line:
                                parts = line.split('\t')
                                if len(parts) > 1 and (parts[1].startswith('ko') or parts[1].startswith('map')):
                                    ko_pathways.append(parts[1])
                                    pathways.add(parts[1])
                    new_ko_entries[ko_id] = ','.join(ko_pathways)
                    sleep(1)  # Respect API rate limits
                except Exception as e:
                    print(f"Error linking KO {k_number} to pathways: {e}")
                    new_ko_entries[ko_id] = ''
        
        for pathway in pathways:
            pathway_genes[pathway].append(row['Gene'])
        gene_pathways[row['Gene']] = list(pathways)
    
    # Save new KO-to-pathway mappings to cache
    if new_ko_entries:
        try:
            new_ko_cache = pd.DataFrame({
                'ko_id': list(new_ko_entries.keys()),
                'pathways': list(new_ko_entries.values())
            })
            if os.path.exists(ko_cache_file):
                existing_ko_cache = pd.read_csv(ko_cache_file)
                new_ko_cache = pd.concat([existing_ko_cache, new_ko_cache]).drop_duplicates('ko_id')
            new_ko_cache.to_csv(ko_cache_file, index=False)
            print(f"Cached {len(new_ko_entries)} new KO-to-pathway mappings to {ko_cache_file}")
        except Exception as e:
            print(f"Error writing to KO cache file: {e}")
    
    print(f"Fetching names for {len(pathway_genes)} pathways")
    pathway_names = get_kegg_pathway_names(pathway_genes.keys(), output_dir)
    
    pathway_data = []
    for pathway, genes in pathway_genes.items():
        pathway_data.append({
            'Pathway_ID': pathway,
            'Pathway_Name': pathway_names.get(pathway, pathway),
            'Gene_Count': len(genes),
            'Genes': ', '.join(genes),
            'Gene_List': genes
        })
    
    pathway_df = pd.DataFrame(pathway_data).sort_values('Gene_Count', ascending=False)
    print(f"Generated pathway DataFrame with {len(pathway_df)} pathways")
    
    top_pathways = pathway_df.head(20)
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(x='Gene_Count', y='Pathway_Name', hue='Pathway_Name', data=top_pathways, 
                     palette='magma', legend=False)
    
    for p in ax.patches:
        width = p.get_width()
        plt.text(width + 0.5, p.get_y() + p.get_height()/2.,
                f'{int(width)}', ha='left', va='center')
    
    plt.title('Top 20 Metabolic Pathways by Gene Count', pad=20)
    plt.xlabel('Number of Genes')
    plt.ylabel('Pathway')
    plt.tight_layout()
    top_pathways_plot_path = os.path.join(output_dir, 'top_pathways.png')
    plt.savefig(top_pathways_plot_path, bbox_inches='tight')
    plt.close()
    print(f"Top pathways plot saved to {top_pathways_plot_path}")
    
    all_pathways = pathway_df['Pathway_ID'].tolist()
    all_genes = list(gene_pathways.keys())
    
    top_pathway_ids = pathway_df.head(15)['Pathway_ID'].tolist()
    matrix_data = []
    
    for gene in all_genes:
        row = {'Gene': gene}
        for pathway in top_pathway_ids:
            row[pathway_names.get(pathway, pathway)] = 1 if pathway in gene_pathways[gene] else 0
        matrix_data.append(row)
    
    matrix_df = pd.DataFrame(matrix_data).set_index('Gene')
    print("Generated gene-pathway matrix")
    
    plt.figure(figsize=(15, 20))
    sns.heatmap(matrix_df.T, cmap='Blues', cbar=False)
    plt.title('Gene-Pathway Association Matrix', pad=20)
    plt.xlabel('Genes')
    plt.ylabel('Pathways')
    plt.xticks([])
    plt.tight_layout()
    heatmap_plot_path = os.path.join(output_dir, 'pathway_heatmap.png')
    plt.savefig(heatmap_plot_path, bbox_inches='tight')
    plt.close()
    print(f"Pathway heatmap saved to {heatmap_plot_path}")
    
    if len(top_pathway_ids) >= 3:
        plt.figure(figsize=(8, 6))
        venn3([set(pathway_genes[top_pathway_ids[0]]),
               set(pathway_genes[top_pathway_ids[1]]),
               set(pathway_genes[top_pathway_ids[2]])],
              (pathway_names.get(top_pathway_ids[0], '...')[:20] + '...',
               pathway_names.get(top_pathway_ids[1], '...')[:20] + '...',
               pathway_names.get(top_pathway_ids[2], '...')[:20] + '...'))
        plt.title('Overlap Between Top 3 Pathways', pad=20)
        overlap_plot_path = os.path.join(output_dir, 'pathway_overlap.png')
        plt.savefig(overlap_plot_path, bbox_inches='tight')
        plt.close()
        print(f"Pathway overlap plot saved to {overlap_plot_path}")
    
    return pathway_df, gene_pathways, matrix_df

def save_results(output_dir, cog_df, pathway_df, gene_pathways, matrix_df):
    """Save all results to Excel files"""
    print("Saving results to Excel...")
    excel_path = os.path.join(output_dir, 'analysis_results.xlsx')
    writer = pd.ExcelWriter(excel_path, engine='xlsxwriter')
    
    cog_df.to_excel(writer, sheet_name='COG_Distribution', index=False)
    pathway_df.to_excel(writer, sheet_name='Pathway_Summary', index=False)
    
    gene_pathway_details = []
    for gene, pathways in gene_pathways.items():
        gene_pathway_details.append({
            'Gene': gene,
            'Pathways': ', '.join(pathways),
            'Pathway_Count': len(pathways)
        })
    pd.DataFrame(gene_pathway_details).to_excel(writer, sheet_name='Gene_Pathway_Mapping', index=False)
    
    matrix_df.reset_index().to_excel(writer, sheet_name='Pathway_Matrix', index=False)
    
    workbook = writer.book
    existing_sheets = set([sheet.name.lower() for sheet in workbook.worksheets()])
    
    images = ['cog_distribution.png', 'top_pathways.png', 'pathway_heatmap.png', 'pathway_overlap.png']
    
    for img in images:
        img_path = os.path.join(output_dir, img)
        if os.path.exists(img_path):
            base_name = img.split('.')[0][:28]
            sheet_name = base_name
            suffix = 1
            while sheet_name.lower() in existing_sheets:
                sheet_name = f"{base_name}_{suffix}"[:31]
                suffix += 1
            try:
                worksheet = workbook.add_worksheet(sheet_name)
                worksheet.insert_image('A1', img_path)
                existing_sheets.add(sheet_name.lower())
            except Exception as e:
                print(f"Error adding worksheet for {img}: {e}")
    
    writer.close()
    print(f"Excel file saved to {excel_path}")

def main():
    parser = argparse.ArgumentParser(description="Analyze functional annotations and generate visualizations.")
    parser.add_argument('input_file', nargs='?', default='sample.annotations.txt', help='Input annotation file (default: sample.annotations.txt)')
    parser.add_argument('--output-dir', default='annotation_analysis_results', help='Output directory for results')
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    print("Parsing annotation file...")
    df = parse_annotations(args.input_file)
    
    print("Analyzing COG categories...")
    cog_df = plot_cog_distribution(df, output_dir)
    
    print("Analyzing metabolic pathways...")
    pathway_df, gene_pathways, matrix_df = analyze_pathways(df, output_dir, use_ko_mapping=True)
    
    print("Saving results...")
    save_results(output_dir, cog_df, pathway_df, gene_pathways, matrix_df)
    
    print(f"\nAnalysis complete! Results saved to: {output_dir}")

if __name__ == "__main__":
    main()