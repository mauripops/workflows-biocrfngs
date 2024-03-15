

rule count_matrix:
    input:
        expand("results/star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "results/counts/all.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="results/counts/all.tsv"
    output:
        "results/deseq2/all.rds"
    params:
        samples=config["samples"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule pca:
    input:
        "results/deseq2/all.rds"
    output:
        report("results/pca.svg", "../report/pca.rst", category="QC", subcategory="PCA")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/pca.log"
    script:
        "../scripts/plot-pca.R"

rule deseq2:
    input:
        "results/deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst", category="Differential Expression", subcategory="{contrast}"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst", category="Differential Expression", subcategory="{contrast}"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"
