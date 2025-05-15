#!/usr/bin/env python3


# Generate barplot of indel rates in Illumina data
# Input: `samtools stats` file
# Output: barplot of indel rates per sample
# Author: Moe


import glob
import sys
import os
import pprint
import subprocess
import io
import matplotlib.pyplot as plt
import pandas as pd


print(sys.executable)
# target_dir = ""
target_dir = sys.argv[1]
files = {os.path.basename(file).replace(".stats.txt", ""): file for file in glob.glob(target_dir + "/*.stats.txt")}
print(f"Input files:")
pprint.pprint(files)

plot_df = {'sample': list(), 'ins_rate': list(), 'del_rate': list()}
for sample, file in files.items():
	results = subprocess.run(f"cat {file} | grep ^SN | cut -f 2-", capture_output=True, text=True, shell=True)
	df = pd.read_csv(io.StringIO(results.stdout), sep="\t", header=None)
	total = float(df[df[0].str.startswith('bases mapped (cigar):')][1].iloc[0])

	results = subprocess.run(f"cat {file} | grep ^ID | cut -f 2-", capture_output=True, text=True, shell=True)
	df = pd.read_csv(io.StringIO(results.stdout), sep="\t", header=None, names=["len", "ins_count", "del_count"])
	df['ins_err'] = df['len'] * df['ins_count'] / total
	df['del_err'] = df['len'] * df['del_count'] / total
	df['sample'] = sample
	df.to_csv(f"{target_dir}/{sample}_indel_rate.tsv", sep="\t", index=None)

	plot_df['sample'].append(sample)
	plot_df['ins_rate'].append(sum(df['ins_err'])) 
	plot_df['del_rate'].append(sum(df['del_err']))

plot_df = pd.DataFrame(plot_df)
plot_df.to_csv(f"{target_dir}/total_indel_rates.tsv", sep="\t", index=None)