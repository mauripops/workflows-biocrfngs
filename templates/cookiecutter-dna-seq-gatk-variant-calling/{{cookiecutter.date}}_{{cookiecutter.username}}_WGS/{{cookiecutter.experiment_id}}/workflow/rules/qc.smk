# Might want to debug
# rule fastqc:
#     input:
#         unpack(get_fastq),
#     output:
#         html="results/qc/fastqc/{sample}-{unit}.html",
#         zip="results/qc/fastqc/{sample}-{unit}.zip",
#     log:
#         "logs/fastqc/{sample}-{unit}.log",
#     wrapper:
#         "v3.4.0/bio/fastqc"

rule samtools_stats:
    input:
        "results/recal/{sample}-{unit}.bam",
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt",
    log:
        "logs/samtools-stats/{sample}-{unit}.log",
    threads: 8
    wrapper:
        "v3.4.0/bio/samtools/stats"


rule multiqc:
    input:
        expand("results/qc/samtools-stats/{u.sample}-{u.unit}.txt", u=units.itertuples())+\
        expand("results/qc/dedup/{u.sample}-{u.unit}.metrics.txt", u=units.itertuples()),
        #expand("results/qc/fastqc/{u.sample}-{u.unit}.zip", u=units.itertuples())+\ TODO debug further once a finished thing is aroung
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.4.0/bio/multiqc"
