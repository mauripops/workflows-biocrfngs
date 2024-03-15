import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from os.path import exists

min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"
#TODO UNCOMMENT & edit SCHEMAS ALL 3
#validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
# validate(units, schema="../schemas/units.schema.yaml")


#### Reference ####
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
config["ref"]["path"] = {}

def valid_paths(paths):
    for path in paths:
        if isinstance(path, str):
            if not exists(path):
                print(f"{path} is not found, check from /data2/biocrfngs/reference_genomes/ or symlink at ./data/")
                return False
        else:
            try:
                iterator = iter(path)
            except TypeError:
                print(path)
                print(f"Path is not a string and not iterable, check what is going on.")
                return False
            else:
                for p in path:
                    if not exists(p):
                        print(f"{p} is not found, check from /data2/biocrfngs/reference_genomes/ or symlink at ./data/")
                        return False
    return True
##### Reference paths #####
config["ref"]["path"]["genome_fasta"] = f"data/{species}/ensembl/{build}/{release}/genome.fasta".replace(' ','')
config["ref"]["path"]["genome_dict"] = f"data/{species}/ensembl/{build}/{release}/genome.dict".replace(' ','')
config["ref"]["path"]["noiupac_variation"] = f"data/{species}/ensembl/{build}/{release}/variation.noiupac.vcf.gz".replace(' ','')
config["ref"]["path"]["tabix_variation"] = f"data/{species}/ensembl/{build}/{release}/genome.fasta".replace(' ','')
config["ref"]["path"]["bwa_index"] = multiext(f"data/{species}/ensembl/{build}/{release}/bwa_mem2/genome.fasta".replace(' ',''), ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
config["ref"]["path"]["vep_cache"] = f"data/{species}/ensembl/{build}/{release}/vep/cache".replace(' ','')
config["ref"]["path"]["vep_plugins"] = f"data/{species}/ensembl/{build}/{release}/vep/plugins".replace(' ','')
config["ref"]["path"]["fai"] = f"data/{species}/ensembl/{build}/{release}/genome.fasta.fai".replace(' ','')
if not valid_paths(config["ref"]["path"].values()):
    print("Reference genomes not all data was found.")
    ## To force error at the moment.
    config = None


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Helper functions #####

# contigs in reference genome
def get_contigs():
    fai = config["ref"]["path"]["fai"]
    #with open(config["ref"]["path"]["fai"],'r') as fai:
    return pd.read_table(fai, header=None, usecols=[0], dtype=str).squeeze()


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}#.values()
    return {"r1": fastqs.fq1}#.values()


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=config["platform"],
    ) +  config["params"]["mapping"]["extra"]


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}-{unit}.fastq".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}-{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def get_regions_param(extra,regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {} ".format(padding)
        return params + extra
    return default + extra


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(extra="",
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        ) + config["params"]["gatk"]["HaplotypeCaller"]["extra"]
    )


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f

#THIS IS NOT USED!
def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL"
    ) + config["params"]["gatk"]["SelectVariants"]["extra"]


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}
