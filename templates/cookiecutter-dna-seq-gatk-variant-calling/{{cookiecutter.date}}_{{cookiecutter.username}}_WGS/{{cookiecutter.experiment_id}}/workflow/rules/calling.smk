if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["path"]["genome_fasta"],
        idx=config["ref"]["path"]["genome_dict"],
        known=config["ref"]["path"]["noiupac_variation"],
        tbi=config["ref"]["path"]["tabix_variation"],
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("results/called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params, #+ config["params"]["gatk"]["HaplotypeCaller"]["extra"],
        java_opts=config["params"]["gatk"]["HaplotypeCaller"]["java"],
    threads: 8#OMP_NUM_THREADS needs to be set as env variable to multithread
    resources:
        mem_mb=32768,
    wrapper:
        "v3.4.0/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["path"]["genome_fasta"],
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf="results/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    params:
        extra=config["params"]["gatk"]["CombineGvcfs"]["extra"],
        java_opts=config["params"]["gatk"]["CombineGvcfs"]["extra"],
    wrapper:
        "v3.4.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["path"]["genome_fasta"],
        gvcf="results/called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]["extra"],
        java_opts=config["params"]["gatk"]["GenotypeGVCFs"]["java"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "v3.4.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    params:
        extra=config["params"]["picard"]["MergeVcfs"]["merge_variants"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    wrapper:
        "v3.4.0/bio/picard/mergevcfs"
