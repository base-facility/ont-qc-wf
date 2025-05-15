#!/usr/bin/env bash


# Generate mpileup output
# Author: Moe


INFILE=$1
OUTFILE=$2

samtools mpileup \
--no-BAQ \
--min-BQ 0 \
--fasta-ref /home/ref/misc/nanoluc_full/phase2_nanoluc_gblock_full_construct.fa \
-o ${OUTFILE} \
--bam-list ${INFILE}