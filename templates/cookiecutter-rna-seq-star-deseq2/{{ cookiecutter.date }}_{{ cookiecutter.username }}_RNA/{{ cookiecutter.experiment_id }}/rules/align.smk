def get_fq(wildcards, read_num=1):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1" if read_num == 1 else "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        return expand("results/trimmed/{sample}-{unit}.{read_num}.fastq.gz",
                      read_num=read_num, **wildcards)

def get_fq1(wildcards):
    return get_fq(wildcards, read_num=1)

def get_fq2(wildcards):
    return get_fq(wildcards, read_num=2)

rule align:
    input:
        fq1 = get_fq1,
    	fq2 = get_fq2,
        idx=config["ref"]["index"],
    output:
        # see STAR manual for additional output files
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}-{unit}/Log.out",
        log_final="results/star/{sample}-{unit}/Log.final.out",
        sj="results/star/{sample}-{unit}/SJ.out.tab",
        unmapped=["results/star/{sample}-{unit}/unmapped.1.fastq.gz","results/star/{sample}-{unit}/unmapped.2.fastq.gz"]
    shadow: "shallow"
    log:
        "results/logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    resources:
        mem_mb=360000
    threads: 40
    wrapper:
        "v3.3.6/bio/star/align"

rule samtools_index:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai"
    params:
        "" # optional params string
    threads:20
    resources:
        mem_mb=180000
    wrapper:
        "v3.3.6/bio/samtools/index"

rule signal_track:
    """
    Create a normalized signal track for display on the UCSC genome browser.
    """
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/tracks/signal/{sample}-{unit}.bw"
    log:
        "results/logs/tracks/signal/{sample}-{unit}.log"
    threads: 20
    resources:
        mem_mb=180000
    conda:
        "../envs/deeptools.yaml"
    shell: "bamCoverage --normalizeUsing RPKM --verbose -b {input[0]} -o {output[0]} -p {threads} &> {log}"
