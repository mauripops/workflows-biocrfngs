rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db")
    log:
        "results/logs/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.Interact.bed"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/qc/rseqc/{sample}-{unit}.junctionanno"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log"
    params:
        extra=r"-q 255",
        prefix="results/qc/rseqc/{sample}-{unit}.junctionsat"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}-{unit}.stats.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_stat/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule rseqc_infer:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.infer_experiment.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_infer/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_innerdis/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.readdistribution.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readdis/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readdup/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readgc/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
    """
    Create a MultiQC report file using the various rule outputs.
    """
    input:
        expand("results/star/{unit.sample}-{unit.unit}/Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        expand("results/logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples())
    output:
        "results/qc/multiqc_report.html",
        directory("results/qc/multiqc_data"),
#        report("results/qc/multiqc_report.html", "../report/multiqc.rst", category="QC", subcategory="MultiQC")
    params:
        extra="--data-dir",
    log:
        "results/logs/multiqc.log"
    params:
        "" # Optional: extra parameters for multiqc.
    wrapper:
        "v3.3.6/bio/multiqc"
