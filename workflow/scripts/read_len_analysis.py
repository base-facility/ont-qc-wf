#!/usr/bin/env python3


# Get read length distribution
# Input: Mapped bam file
# Output: Read length distribution histogram plots

import pysam
import matplotlib.pyplot as plt
import numpy as np
import sys
from bioUtils import seqUtils

infile = sys.argv[1]
outprefix = sys.argv[2]
sample = sys.argv[3]

def plot_read_length_histogram(bam_file, output_prefix, sample_name):
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        # read_lengths = [read.query_length for read in bam if read.is_mapped]
        read_length = []
        for read in bam:
            if read.is_mapped  and not read.is_supplementary and not read.is_duplicate and not read.is_secondary and read.mapping_quality == 60:
                query = read.query_sequence
                query_id = read.query_name
                # Debug
                if not query:
                    print(query_id)
                # Without removing polyA
                # read_length.append(len(query))
                # With polyA removed
                seq_obj = seqUtils(query)
                polished_query = seq_obj.rmPolyA()
                read_length.append(len(polished_query))


    # Plot histogram of read lengths
    plt.figure(figsize=(10, 6))
    xlim = 2000
    bins = np.arange(0, xlim + 10, 10)
    plt.hist(read_length, bins=bins, color='skyblue', edgecolor='black')
    plt.xlim(0, xlim)
    plt.title(f'{sample_name} - Distribution of Read Lengths')
    plt.xlabel('Read Length (bp)')
    plt.ylabel('Frequency')
    
    # Save the plot in different formats
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.svg", format='svg')  # Vector format (SVG)
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.png", format='png', dpi=300)  # Rasterized format (PNG)
    plt.savefig(f"{output_prefix}/{sample_name}_rl_hist.pdf", format='pdf')  # Vector format (PDF)

plot_read_length_histogram(infile, outprefix, sample)

# Dev
# read.get_tag("qs") >= 10
# read.mapping_quality == 60
# seq_obj = seqUtils(query)
# polished_query = seq_obj.rmPolyA()
# read_length.append(len(polished_query))
