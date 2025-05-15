#!/usr/bin/env python3


# Get polya tail length distribution
# Input: Mapped bam file
# Output: box plots of polya tail length distribution

# TODO
# Modify script to make it compatible with Snakemake

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import sys

files = {
    "bicycle_one": "/home/uqmzardb/net/rmt_gridion/projects/bicycle/250416_ont_drna_data/base_calls/bicycle_one/calls_2025-04-21_T23-41-45.bam",
    "bicycle_sep": "/home/uqmzardb/net/rmt_gridion/projects/bicycle/250416_ont_drna_data/base_calls/bicycle_sep/calls_2025-04-17_T06-02-06.bam",
    "pcr": "/home/uqmzardb/net/rmt_gridion/projects/bicycle/250416_ont_drna_data/base_calls/pcr/calls_2025-04-17_T07-14-35.bam",
    "plasmid": "/home/uqmzardb/net/rmt_gridion/projects/bicycle/250416_ont_drna_data/base_calls/plasmid/calls_2025-04-17_T09-05-44.bam",
    # "refstd_126a": "/home/uqmzardb/net/rmt_gridion/projects/bicycle/250416_ont_drna_data/base_calls/refstd_126a/calls_2025-04-17_T10-44-58.bam",
}
output_prefix = "/home/uqmzardb/projects/bicycle-qc-wf/results/polya_len_plot"

def get_pt_df(infile):
    with pysam.AlignmentFile(infile, 'rb', check_sq=False) as bam:
    # Total number of reads
        total_reads = 0
        records = {'id': [],
                'pt': []}
        for read in bam:
            # Only primary mapped reads with pt tag and high mapping quality
            if read.is_mapped  and not read.is_supplementary and not read.is_duplicate and not read.is_secondary and read.mapping_quality == 60:
                total_reads += 1
                if read.has_tag('pt'):
                    records['id'].append(read.query_name)
                    records['pt'].append(int(read.get_tag('pt')))
        df = pd.DataFrame(records)
        print(f"polya ratio: {len(df['pt']) / total_reads}")
        return((df, total_reads))

output_df = pd.DataFrame({'id': [], 
                          'pt': [],
                          'sample': []})
ratio_df = {'sample': [],
            'polya_ratio': []}
for sample, infile in files.items():
    print(f"sample: {sample}")
    _out = get_pt_df(infile)
    df = _out[0]
    total_reads = _out[1]
    # Sanity check
    print(f"{sample} total reads: {total_reads}")
    ratio_df['sample'].append(sample)
    ratio_df['polya_ratio'].append(len(df['pt']) / total_reads)
    df['sample'] = sample
    output_df = pd.concat([output_df, df], ignore_index=True)

# Prepare data for boxplot
sample_names = output_df['sample'].unique()
data = [output_df[output_df['sample'] == sample]['pt'] for sample in sample_names]

# Create boxplot
plt.figure(figsize=(10, 7))
# plt.axhline(y=100, color='b', linestyle='--', linewidth=1)
plt.axhline(y=126, color='b', linestyle='--', linewidth=1)
plt.boxplot(data, labels=sample_names, showfliers=False)
plt.title('Distribution of polyA/T Tail Lengths')
plt.xlabel('Sample')
plt.ylabel('Length (bp)')
plt.xticks(rotation=45)
# plt.grid(axis='y', linestyle='--', alpha=0.7)

plt.savefig(f"{output_prefix}/polya_len_boxplot.svg", format="svg")
plt.savefig(f"{output_prefix}/polya_len_boxplot.png", dpi=300)
plt.savefig(f"{output_prefix}/polya_len_boxplot.pdf", format='pdf')

# Geneerate barplot
fig, ax = plt.subplots(figsize=(8, 6))
bars = ax.bar(ratio_df['sample'], ratio_df['polya_ratio'], 
              color=['blue', 'green', 'red'])

# Add text labels on top of each bar
for bar in bars:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom')

# Add labels and title
ax.set_xlabel('Sample')
ax.set_ylabel('Polya Ratio')
ax.set_title('PolyA ratios per sample')

fig.savefig(f"{output_prefix}/polya_ratios_barplot.svg", format="svg")
fig.savefig(f"{output_prefix}/polya_ratios_barplot.png", dpi=300)
fig.savefig(f"{output_prefix}/polya_ratios_barplot.pdf", format='pdf')

# Dev
# read.mapping_quality == 60
# read.get_tag("qs") >= 8