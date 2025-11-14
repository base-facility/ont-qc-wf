#!/usr/bin/env python3


# Get read length distribution
# Input: Mapped bam file
# Output: Read length distribution histogram plots
# Author: Moe

import pysam, json
import matplotlib.pyplot as plt
import numpy as np
import sys
from bioUtils import seqUtils
import statistics

infile = sys.argv[1]
outprefix = sys.argv[2]
sample = sys.argv[3]

def plot_read_length_histogram(bam_file, output_prefix, sample_name):
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        discarded_ids = []
        alignment_spans = []
        read_length = []
        # filter primary mapped reads only
        for read in bam:
            if (not read.query_sequence or read.is_supplementary or read.is_duplicate or read.is_secondary or read.is_qcfail or read.is_unmapped):
                discarded_ids.append(read.query_name)
                continue
            if read.is_mapped:
                query = read.query_sequence

                # alignment span as length
                span_len = read.reference_end - read.reference_start + 1
                alignment_spans.append(span_len)
                
                # Without removing polyA
                # read_length.append(len(query))
                
                # With polyA removed
                seq_obj = seqUtils(query)
                polished_query = seq_obj.rmPolyA()
                read_length.append(len(polished_query))
            else:
                print("*** ERROR ***")
                print(read.query_name)
                break

        json.dump(discarded_ids, open(f'{output_prefix}/{sample_name}_discarded_ids.json', 'w'), indent=4)
        json.dump(alignment_spans, open(f'{output_prefix}/{sample_name}_align_spans.json', 'w'), indent=4)
        json.dump(read_length, open(f'{output_prefix}/{sample_name}_rl.json', 'w'), indent=4)


    # Plot histogram of read lengths
    plt.figure(figsize=(10, 6))
    # set xlim to 99th percentile to remove very long concatermer from plot
    xlim = np.percentile(read_length, 99)
    # clipped = np.clip(read_length, None, xlim)
    bins = np.arange(0, xlim + 10, 10)
    median_rl = np.median(read_length)
    mean_rl = np.mean(read_length)
    plt.axvline(median_rl, color='red', linestyle='--', label=f"Median: {median_rl:.0f}")
    plt.axvline(mean_rl, color='green', linestyle='-', label=f"Mean: {mean_rl:.0f}")

    plt.hist(read_length, bins=bins, color='skyblue', edgecolor='black')
    plt.xlim(0, xlim + 10)
    plt.title(f'{sample_name} - Distribution of Read Lengths')
    plt.xlabel('Read Length (bp)')
    plt.ylabel('Frequency')
    
    # Save the plot in different formats
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.svg", format='svg')  # Vector format (SVG)
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.png", format='png', dpi=300)  # Rasterized format (PNG)
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.pdf", format='pdf')  # Vector format (PDF)

    plt.close()
    
    # Plot histogram of read spans
    plt.figure(figsize=(10, 6))
    bins = np.arange(0, max(alignment_spans) + 10, 10)
    plt.hist(alignment_spans, bins=bins, color='skyblue', edgecolor='black')
    plt.title(f'{sample_name} - Distribution of alignment lengths')
    plt.xlabel('Read length (bp)')
    plt.ylabel('Frequency')
    
    # Save the plot in different formats
    plt.savefig(f"{output_prefix}/{sample_name}_alignment_span_hist.svg", format='svg')  # Vector format (SVG)
    plt.savefig(f"{output_prefix}/{sample_name}_alignment_span_hist.png", format='png', dpi=300)  # Rasterized format (PNG)
    plt.savefig(f"{output_prefix}/{sample_name}_alignment_span_hist.pdf", format='pdf')  # Vector format (PDF)

    plt.close()

plot_read_length_histogram(infile, outprefix, sample)