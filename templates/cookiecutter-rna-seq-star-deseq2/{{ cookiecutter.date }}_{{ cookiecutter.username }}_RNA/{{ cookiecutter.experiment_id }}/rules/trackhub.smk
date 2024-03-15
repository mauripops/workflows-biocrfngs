localrules: get_genome_sizes, trackhub

rule get_genome_sizes:
    output:
        temp("results/chrom.sizes")
    params:
        genome=config["ref"]["version"]
    conda:
        "../envs/ucsc-bedtobigbed.yaml"
    shell:
        "export LC_COLLATE=C; "
        "fetchChromSizes {params.genome} | sort -k1,1 -k2,2n > {output[0]}"

rule fix_interact_scores:
    """
    Scores in BED need to be between 0-1000 or we can't convert them to bigBed.
    """
    input:
        bed="results/qc/rseqc/{sample}-{unit}.junctionanno.junction.Interact.bed",
        sizes="results/chrom.sizes"
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.Interact.fixed_scores.bed"
    wrapper:
        "file:wrappers/fix_peak_scores"

rule interact_to_bigbed:
    """
    Convert the interaction BED file created by RSeQC into a bigBed file for display on the genome browser.
    """
    input:
        bed="results/qc/rseqc/{sample}-{unit}.junctionanno.junction.Interact.fixed_scores.bed",
        sizes="results/chrom.sizes"
    output:
        "results/tracks/interaction/{sample}-{unit}.bb"
    params:
        as_file="scripts/interact.as",
    conda:
        "../envs/ucsc-bedtobigbed.yaml"
    shell:
        "bedToBigBed -as={params.as_file} -type=bed5+13 {input.bed} {input.sizes} {output[0]}"

rule trackhub:
    """
    Automatically create all the files needed for a UCSC track hub.
    """
    input:
        signals=expand("results/tracks/signal/{unit.sample}-{unit.unit}.bw", unit=units.itertuples()),
        interacts=expand("results/tracks/interaction/{unit.sample}-{unit.unit}.bb", unit=units.itertuples())
    output:
        directory("results/trackhub"),
        genomes="results/trackhub/{hub_name}.genomes.txt".format(hub_name=config["trackhub"]["name"]),
        hub="results/trackhub/{hub_name}.hub.txt".format(hub_name=config["trackhub"]["name"])
    log:
        "results/logs/trackhub.log"
    params:
        samples=samples,
        units=units
    conda:
        "../envs/trackhub.yaml"
    script:
        "../scripts/build_trackhub.py"
