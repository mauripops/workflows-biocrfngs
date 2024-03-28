def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

ruleorder: cutadapt_pe > cutadapt

rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
    threads: 40
    params:
        adapters="-a {0} -A {0} {1}".format(config["trimming"]["adapter"], config["params"]["cutadapt-pe"]),
        extra="--minimum-length=30 --nextseq-trim=30 -a A\{100\} -A A\{100\}"
    resources:
        mem_mb=360000
	# This is so that it only keeps a set of reads if both pairs are at least 20bp long
	# And, it performs 3' poly-A trimming on the reads
    log: "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "v3.3.6/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastq
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
    threads: 40
    params:
        adapters="-a {0} {1}".format(config["trimming"]["adapter"], config["params"]["cutadapt-se"]),
        extra="--minimum-length=35 --nextseq-trim=30 -a A\{100\}"
    resources:
        mem_mb=360000
    log: "results/logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "v3.3.6/bio/cutadapt/se"
