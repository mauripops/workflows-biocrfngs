localrules: download_files_local

rule download_files_local:
    """
    Unzip the fastq files from the local server to the project folder.
    """
    output:
        temp("data/reads/{sequencing_run}_{lane}_{barcode}_{read}.fq.gz"),
    shell:
        """
        cp /data3/ztron/autorunDW/DNBSEQ-T7/R*/*/ztron/{wildcards.sequencing_run}_{wildcards.lane}_[0-9][0-9][0-9][0-9][0-9]/{wildcards.sequencing_run}_{wildcards.lane}/`basename {output[0]}` {output[0]}
        """