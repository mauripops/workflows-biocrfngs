.. raw:: html

   <div class="row justify-content-right" id="brand_big_logo"></div>

This report describes the analysis steps and results obtained from running a custom {{ snakemake.config['workflow']['name'] }} Snakemake analysis workflow. This workflow performs differential expression analysis on paired-end RNA-seq data. These results were produced with v{{ snakemake.config['workflow']['version']['number'] }} of the workflow released on {{ snakemake.config['workflow']['version']['date'] }}.

**Analysis Steps**

After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_, reads were mapped to the {{ snakemake.config['ref']['organism_label'] }} {{ snakemake.config['ref']['provider'] }} {{ snakemake.config['ref']['version'] }} reference genome followed by gene count quantification using `STAR <https://github.com/alexdobin/STAR>`_. The gene counts of replicates were summed up. Integrated normalization and differential expression analysis was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ following the standard procedure as outlined in this `tutorial <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html>`_.
