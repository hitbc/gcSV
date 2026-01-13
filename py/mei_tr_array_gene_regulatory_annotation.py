
import gzip
import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# Parse MEI interval file
def parse_mei_regions(file_path):
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Format: chrX_309000_311500
            parts = line.split('_')
            chrom = parts[0]  # Keep 'chr' prefix to match with GTF file
            start = int(parts[1])
            end = int(parts[2])
            regions.append((chrom, start, end))
    return regions

# Parse GFF3 file
def parse_gff3(file_path):
    regulatory_features = []
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            chrom = parts[0]
            # Add chr prefix to chromosome names in GFF3 file to match with MEI regions
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            full_line = line
            regulatory_features.append((chrom, start, end, feature_type, full_line))
    return regulatory_features

# Parse GTF file
def parse_gtf(file_path):
    gene_features = []
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            chrom = parts[0]
            # Chromosome names in GTF file already have chr prefix, no need to process
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            full_line = line
            gene_features.append((chrom, start, end, feature_type, full_line))
    return gene_features

# Check if two intervals overlap
def intervals_overlap(start1, end1, start2, end2):
    return start1 <= end2 and start2 <= end1

# Analyze overlap between MEI intervals and annotation features
def analyze_overlap(mei_regions, all_features):
    # Group annotation features by chromosome
    features_by_chrom = defaultdict(list)
    for feature in all_features:
        chrom, start, end, feature_type,full_line = feature
        features_by_chrom[chrom].append((start, end, feature_type,full_line))

    # Statistical results
    overlap_stats = defaultdict(int)
    total_mei = 0
    overlap_details = []  # Store detailed information for each overlap
    overlapping_mei = set()  # Record MEI intervals that overlap with any feature

    for mei_chrom, mei_start, mei_end in mei_regions:
        total_mei += 1
        if mei_chrom not in features_by_chrom:
            continue

        # Check overlap with all annotation features on this chromosome
        has_overlap = False
        # Record feature types that this MEI interval has overlapped with, ensuring each type is counted only once
        overlapped_types = set()

        for reg_start, reg_end, feature_type, full_line in features_by_chrom[mei_chrom]:
            if intervals_overlap(mei_start, mei_end, reg_start, reg_end):
                # Only count if this MEI interval has not yet overlapped with this type of feature
                if feature_type not in overlapped_types:
                    overlap_stats[feature_type] += 1
                    overlapped_types.add(feature_type)

                has_overlap = True
                # Add detailed overlap information
                overlap_details.append({
                    'mei_region': f"{mei_chrom}_{mei_start}_{mei_end}",
                    'feature_type': feature_type,
                    'feature_region': f"{mei_chrom}_{reg_start}_{reg_end}",
                    'full_line' : full_line
                })

        if has_overlap:
            overlapping_mei.add(f"{mei_chrom}_{mei_start}_{mei_end}")

    # Calculate the number of non-overlapping MEIs
    non_overlapping = total_mei - len(overlapping_mei)

    return total_mei, overlap_stats, overlap_details, non_overlapping

# Save detailed overlap information to file
def save_overlap_details(overlap_details, output_file):
    with open(output_file, 'w') as f:
        # Write header row
        f.write("MEI_Region\tFeature_Type\tFeature_Region\n")

        # Write each overlap record
        for detail in overlap_details:
            f.write(f"{detail['mei_region']}\t{detail['feature_type']}\t{detail['feature_region']}\t{detail['full_line']}\n")

    print(f"Detailed overlap information saved to: {output_file}")

# Generate result table
def generate_table(total_mei, overlap_stats, non_overlapping, output_file):
    # Create DataFrame
    data = []
    for feature_type, count in sorted(overlap_stats.items()):
        percentage = (count / total_mei) * 100
        data.append({
            "Annotation Feature": feature_type,
            "Count": count,
            "Percentage (%)": f"{percentage:.2f}"
        })

    # Add number of MEIs that do not fall on annotation features
    non_overlapping_percentage = (non_overlapping / total_mei) * 100
    data.append({
        "Annotation Feature": "No overlap",
        "Count": non_overlapping,
        "Percentage (%)": f"{non_overlapping_percentage:.2f}"
    })

    df = pd.DataFrame(data)

    # Save as CSV file
    df.to_csv(output_file, index=False)
    print(f"Result table saved to: {output_file}")

    return df

# Generate bar chart
def generate_bar_chart(total_mei, overlap_stats, non_overlapping, output_file):
    # Define the order of feature types
    feature_order = ['gene', 'exon', 'CDS', 'UTR', 
                    'enhancer', 'promoter', 'CTCF_binding_site', 'open_chromatin_region']

    # Prepare data
    feature_types = []
    counts = []

    # Add feature types and counts in the specified order
    for feature_type in feature_order:
        if feature_type == 'No overlap':
            counts.append(non_overlapping)
        elif feature_type in overlap_stats:
            counts.append(overlap_stats[feature_type])
        else:
            continue  # Skip if feature type does not exist
        feature_types.append(feature_type)

    # Create bar chart
    plt.figure(figsize=(12, 8))
    # Use more colors to accommodate potentially more feature types
    colors = plt.cm.tab20.colors
    bars = plt.bar(feature_types, counts, color=colors[:len(feature_types)])

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{int(height)}',
                 ha='center', va='bottom')

    # Set chart title and labels
    plt.title('MEI Intervals Overlapping with Annotation Features', fontsize=14)
    plt.xlabel('Annotation Feature Type', fontsize=12)
    plt.ylabel('Number of MEI Intervals', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save chart
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Bar chart saved to: {output_file}")
    plt.close()

def process_mei_file(mei_file):
    # Get the directory and base name of the input file
    input_dir = os.path.dirname(mei_file)
    base_name = os.path.basename(mei_file)
    name_without_ext = os.path.splitext(base_name)[0]
    
    # Define output file paths, using input file path + unique suffix
    table_output = os.path.join(input_dir, f"{name_without_ext}_annotation_overlap.csv")
    chart_output = os.path.join(input_dir, f"{name_without_ext}_annotation_overlap.png")
    details_output = os.path.join(input_dir, f"{name_without_ext}_annotation_details.tsv")
    
    # Now these output paths can be used
    print("Input file:", mei_file)
    print("Table output:", table_output)
    print("Chart output:", chart_output)
    print("Details output:", details_output)
    
    return table_output, chart_output, details_output

import os
import argparse

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Process MEI file and generate output files')
    
    # Add input file parameter
    parser.add_argument('mei_file', type=str, help='Path to the input MEI file')
    
    # Parse command line arguments
    args = parser.parse_args()
    
    # Process input file
    table_output, chart_output, details_output = process_mei_file(args.mei_file)
    
    gff3_file = 'Homo_sapiens.GRCh38.regulatory_features.v115.gff3.gz'
    gtf_file = 'gencode.v44.basic.annotation.sorted.gtf.gz'

    # Parse files
    print("Parsing MEI interval file...")
    mei_regions = parse_mei_regions(args.mei_file)

    print("Parsing GFF3 annotation file...")
    regulatory_features = parse_gff3(gff3_file)

    print("Parsing GTF annotation file...")
    gene_features = parse_gtf(gtf_file)

    # Merge all annotation features
    all_features = regulatory_features + gene_features

    # Analyze overlap
    print("Analyzing overlap between MEI intervals and all annotation features...")
    total_mei, overlap_stats, overlap_details, non_overlapping = analyze_overlap(mei_regions, all_features)

    # Save detailed overlap information
    save_overlap_details(overlap_details, details_output)

    # Output results
    print(f"\nTotal analyzed {total_mei} MEI intervals")
    print("\nNumber of MEIs falling on different types of annotation features:")
    for feature_type, count in sorted(overlap_stats.items()):
        print(f"{feature_type}: {count}")

    # Calculate proportions
    print("\nProportion of each type of annotation feature:")
    for feature_type, count in sorted(overlap_stats.items()):
        percentage = (count / total_mei) * 100
        print(f"{feature_type}: {percentage:.2f}%")

    # Generate table and chart
    print("\nGenerating result table...")
    df = generate_table(total_mei, overlap_stats, non_overlapping, table_output)
    print(df)

    print("\nGenerating bar chart...")
    generate_bar_chart(total_mei, overlap_stats, non_overlapping, chart_output)

if __name__ == "__main__":
    main()