| `MA plot <https://en.wikipedia.org/wiki/MA_plot>`_ of log fold change vs. mean of normalized counts for each gene when calculating the differential expression for contrast {{ snakemake.wildcards.contrast }}.
|
| The log2FoldChange values represent log2 ( {{ snakemake.params.contrast[0] }} / {{ snakemake.params.contrast[1] }} ).
|
| Please see `this section <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#ma-plot>`_ of the DESeq2 documentation for more information.
