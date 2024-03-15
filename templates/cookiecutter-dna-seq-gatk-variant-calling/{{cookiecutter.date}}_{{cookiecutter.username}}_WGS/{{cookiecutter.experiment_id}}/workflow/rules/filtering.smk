rule select_calls:
    input:
        ref=config["ref"]["path"]["genome_fasta"],
        vcf="results/genotyped/all.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.vcf.gz"),
    params:
        extra=get_vartype_arg, #+ config["params"]["gatk"]["SelectVariants"]["extra"],
        java_opts=config["params"]["gatk"]["SelectVariants"]["java"],
    log:
        "logs/gatk/selectvariants/{vartype}.log",
    resources:
        mem_mb=32768,
    threads: 4
    wrapper:
        "v3.4.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref=config["ref"]["path"]["genome_fasta"],
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.hardfiltered.vcf.gz"),
    params:
        filters=get_filter,
        extra=config["params"]["gatk"]["HardFiltering"]["extra"],
        java_opts=config["params"]["gatk"]["HardFiltering"]["java"],
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    resources:
        mem_mb=32768,
    threads: 4
    wrapper:
        "v3.4.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.recalibrated.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]["extra"],
        java_opts=config["params"]["gatk"]["VariantRecalibrator"]["java"],
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log",
    resources:
        mem_mb=32768,
    threads: 4
    wrapper:
        "v3.4.0/bio/gatk/variantrecalibrator"


rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    params:
        extra=config["params"]["picard"]["MergeVcfs"]["merge_calls"],
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=32768,
    threads: 8
    wrapper:
        "v3.4.0/bio/picard/mergevcfs"
