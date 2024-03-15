localrules: download_files_local

rule download_files_local:
    """
    Copy the fastq files from the local server to the project folder.
    """
    output:
        protected("resources/reads/{sequencing_run}_{lane}_{barcode}_1.fq.gz"),
        protected("resources/reads/{sequencing_run}_{lane}_{barcode}_2.fq.gz")
    wildcard_constraints:
        lane="L01|L02|L03|L04",
        barcode="\d+"
    shell:
        "cp /disk2/bgi/R1100400190016B/{wildcards.sequencing_run}/{wildcards.lane}/`basename {output[0]}` {output[0]}; "
        "cp /disk2/bgi/R1100400190016B/{wildcards.sequencing_run}/{wildcards.lane}/`basename {output[1]}` {output[1]}; "

#rule download_files:
#    """
#    Download the fastq files from the remote server using rsync.
#    """
#    output:
#        protected("resources/reads/{sequencing_run}_{lane}_{barcode}_1.fq.gz"),
#        protected("resources/reads/{sequencing_run}_{lane}_{barcode}_2.fq.gz")
#    wildcard_constraints:
#        lane="L01|L02|L03|L04",
#        barcode="\d+"
#    shell:
#        "rsync -v -h --progress --archive --no-D ngs.ust.hk:/disk2/bgi/R1100400190016B/{wildcards.sequencing_run}/{wildcards.lane}/`basename {output[0]}` {output[0]}; "
#        "rsync -v -h --progress --archive --no-D ngs.ust.hk:/disk2/bgi/R1100400190016B/{wildcards.sequencing_run}/{wildcards.lane}/`basename {output[1]}` {output[1]}; "
