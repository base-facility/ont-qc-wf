# Convert bam to fastq
if not config['samples_mapped'] and not config['fastq_input']:
    rule bam2fastq:
        input: lambda wildcards: input_files[wildcards.id]
        output: "results/fastq/{id}.fastq"
        conda: "../envs/samtools.yaml"
        log: "logs/{id}_bam2fastq.log"
        benchmark: "benchmarks/{id}_bam2fastq.benchmark"
        threads: 8
        shell:
            '''
            samtools fastq -@ {threads} {input} > {output} 2> {log}
            '''

if not config['samples_mapped']:
# Subsample input fastq files
    rule subsample:
        input: lambda wildcards: input_files[wildcards.id] if config['fastq_input'] else "results/fastq/{id}.fastq"
        output: "results/subsample/{id}_sub.fastq.gz"
        threads: 8
        conda: "../envs/seqtk.yaml"
        log: "logs/subsample_{id}.log"
        benchmark: "benchmarks/subsample_{id}.benchmark"
        params:
            size = 1000000
        shell:
        '''
        seqtk sample -s60 {input} {params.size} | bgzip -@ {threads} > {output}
        '''

if not config['samples_mapped']:
    rule minimap2:
        input: "results/subsample/{id}_sub.fastq.gz"
        output: "results/minimap2/{id}_mapped.sorted.bam",
        conda: "../envs/minimap2.yaml"
        log: "logs/minimap2_{id}.log"
        benchmark: "benchmarks/minimap2_{id}.benchmark"
        threads: 8
        params:
            ref_fa = input_refs[wildcards.id]
        shell:
            '''
            minimap2 -t {threads} \
            -ax map-ont \
            --secondary=no \
            -k 14 \
            {params.ref_fa} \
            {input} | \
            samtools sort -@ {threads} -O BAM -o {output} -;
            samtools index -@ {threads} {output} 2> {log}
            '''

# rule samtools_stats:
#     input: 
#         lambda wildcards: "results/minimap2/{id}_mapped.sorted.bam" if config['fastq_input'] else input_files[wildcards.id]
#     output: "results/samtools_stats/{id}.stats.txt",
#     conda: "../envs/samtools.yaml"
#     log: "logs/samtools_stats_{id}.log"
#     benchmark: "benchmarks/samtools_stats_{id}.benchmark"
#     threads: 8
#     params:
#         prefix = "results/samtools_stats/"
#     shell:
#         '''
#         samtools stats -@ {threads} {input} > {output} 2> {log}
#         '''

# rule mpileup:
#     input: "results/minimap2/{id}_mapped.sorted.bam"
#     output: "results/mpileup/{id}.mpileup.tsv",
#     conda: "../envs/samtools.yaml"
#     log: "logs/samtools_{id}.log"
#     benchmark: "benchmarks/samtools_{id}.benchmark"
#     params:
#         ref_fa = input_refs[wildcards.id]
#     shell:
#         '''
#         samtools mpileup \
#         --no-BAQ \
#         --min-BQ 0 \
#         --fasta-ref {params.ref_fa} \
#         -o {output} \
#         {input} 2> {log}
#         '''

# rule indel_rates:
#     input: "results/samtools_stats/{id}.stats.txt"
#     output: "results/samtools_stats/{id}_indel_rate.tsv",
#     conda: "../envs/python.yaml"
#     log: "logs/indel_rates_{id}.log"
#     benchmark: "benchmarks/indel_rates_{id}.benchmark"
#     params:
#         prefix = "results/samtools_stats/"
#     shell:
#         '''
#         python workflow/scripts/getIndelRate.py {params.prefix} 2> {log}
#         '''

rule samtools_filter:
    input:
        lambda wildcards: "results/minimap2/{id}_mapped.sorted.bam" if config['fastq_input'] else input_files[wildcards.id]
    output: "results/primary_mapped/{id}.primary_mapped.sorted.bam"
    conda: "../envs/samtools.yaml"
    log: "logs/samtools_filter_primary_mapped_{id}.log"
    benchmark: "benchmark/samtools_filter_primary_mapped_{id}.benchmark"
    shell:
        '''
        samtools view -o {output} -O BAM -h -F 0xF04 {input} 2> {log};
        samtools index -@ 8 {input}
        '''

rule samtools_stats:
    input: "results/primary_mapped/{id}.primary_mapped.sorted.bam"
    output: "results/samtools_stats_primary_mapped/{id}.stats.txt",
    conda: "../envs/samtools.yaml"
    log: "logs/samtools_stats_primary_mapped_{id}.log"
    benchmark: "benchmarks/samtools_stats_primary_mapped_{id}.benchmark"
    threads: 8
    shell:
        '''
        samtools stats -@ {threads} {input} > {output} 2> {log}
        '''

rule mpileup:
    input: "results/primary_mapped/{id}.primary_mapped.sorted.bam"
    output: "results/mpileup_primary_mapped/{id}.mpileup.tsv",
    conda: "../envs/samtools.yaml"
    log: "logs/samtools_mpileup_primary_mapped{id}.log"
    benchmark: "benchmarks/samtools_mpileup_primary_mapped_{id}.benchmark"
    params:
        ref_fa = input_refs[wildcards.id]
    shell:
        '''
        samtools mpileup \
        --no-BAQ \
        --min-BQ 0 \
        --fasta-ref {params.ref_fa} \
        -o {output} \
        {input} 2> {log}
        '''

# compute indel rates for primary mapped reads
rule indel_rates:
    input: "results/samtools_stats_primary_mapped/{id}.stats.txt"
    output: "results/samtools_stats_primary_mapped /{id}_indel_rate.tsv",
    conda: "../envs/python.yaml"
    log: "logs/indel_rates_primary_mapped_{id}.log"
    benchmark: "benchmarks/indel_rates_primary_mapped_{id}.benchmark"
    params:
        prefix = "results/samtools_stats_primary_mapped/"
    shell:
        '''
        python workflow/scripts/getIndelRate.py {params.prefix} 2> {log}
        '''

rule read_len_hist:
    input: "results/primary_mapped/{id}.primary_mapped.sorted.bam"
    output: "results/read_len_hist/{id}_rl_hist.svg",
    conda: "../envs/pysam.yaml"
    log: "logs/read_len_hist_{id}.log"
    benchmark: "benchmarks/read_len_hist_{id}.benchmark"
    params:
        prefix = "results/read_len_hist/"
    shell:
        '''
        python workflow/scripts/read_len_analysis.py {input} {params.prefix} {wildcards.id} 2> {log}
        '''