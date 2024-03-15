import os

localrules: to_server

rule to_server:
    output: touch("xfer_upload.done")
    params:
        parent_folder=os.path.basename(os.path.dirname(os.path.realpath(os.getcwd()))),
        workflow_folder=os.path.basename(os.getcwd())
    shell:
        """
        rsync --rsync-path="mkdir -p /disk2/workflows/{params.parent_folder}/{params.workflow_folder} && rsync" -v -h --progress --archive --no-D results resources units.tsv samples.tsv ngs.ust.hk:/disk2/workflows/{params.parent_folder}/{params.workflow_folder}
        """
