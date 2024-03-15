rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        adapters="-a {0} -g {1}".format(config["params"]["trimming"]["adapter-a"], config["params"]["trimming"]["adapter-g"]),
        extra=config["params"]["trimming"]["extra-se"],
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    threads: 2
    wrapper:
        "v3.4.0/bio/cutadapt/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    #Chris: This is so that it only keeps a set of reads if both pairs are at least 20bp long
	# And, it performs 3' poly-A trimming on the reads
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        adapters="-a {0} -A {0} -g {1} -G {1}".format(config["params"]["trimming"]["adapter-a"], config["params"]["trimming"]["adapter-g"]),
        extra=config["params"]["trimming"]["extra-pe"],
    resources:
        mem_mb=72000
    log: "logs/cutadapt/{sample}-{unit}.log"
    threads: 8
    wrapper:
        "v3.4.0/bio/cutadapt/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=config["ref"]["path"]["bwa_index"],
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group, # +  config["params"]["mapping"]["extra"]
        sort="samtools",#var called sorting in wrapper v3.4.0/bio/bwa/mem"
        sort_order="coordinate",
    threads: 20
    wrapper:
        "v3.4.0/bio/bwa-mem2/mem"


rule mark_duplicates:
    input:
        bams="results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=temp("results/dedup/{sample}-{unit}.bam"),
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        extra=config["params"]["picard"]["MarkDuplicates"]["extra"],
        java_opts=config["params"]["picard"]["MarkDuplicates"]["java"],# using 4 GCthreads reason https://www.biorxiv.org/content/10.1101/348565v1.full.pdf 
    resources:
        mem_mb=32768,
    threads: 8 #to prevent GC competing for threads
    wrapper:
        "v3.4.0/bio/picard/markduplicates"

# rule recalibrate_base_qualities:
#     input:
#         bam=get_recal_input(),
#         bai=get_recal_input(bai=True),
#         ref=config["ref"]["path"]["genome_fasta"],
#         dict=config["ref"]["path"]["genome_dict"],
#         known=config["ref"]["path"]["noiupac_variation"],
#         known_idx=config["ref"]["path"]["tabix_variation"],
#     output:
#         recal_table="results/recal/{sample}-{unit}.grp",
#     log:
#         "logs/gatk/bqsr/{sample}-{unit}.log",
#     params:
#         extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
#     resources:
#         mem_mb=1024,
#     wrapper:
#         "v3.4.0/bio/gatk/baserecalibrator"
rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["path"]["genome_fasta"],
        dict=config["ref"]["path"]["genome_dict"],
        known=config["ref"]["path"]["noiupac_variation"],
        known_idx=config["ref"]["path"]["tabix_variation"],
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(config["params"]["gatk"]["BaseRecalibrator"]["extra"]),
        java_opts=config["params"]["gatk"]["BaseRecalibrator"]["java"],
    resources:
        mem_mb=32768,
    threads: 10
    wrapper:
        "v3.4.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["path"]["genome_fasta"],
        dict=config["ref"]["path"]["genome_dict"],
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam=protected("results/recal/{sample}-{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(config["params"]["gatk"]["ApplyBaseQualityRecalibrator"]["extra"]),
        java_opts=config["params"]["gatk"]["ApplyBaseQualityRecalibrator"]["java"],
    wrapper:
        "v3.4.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    threads: 16
    wrapper:
        "v3.4.0/bio/samtools/index"
