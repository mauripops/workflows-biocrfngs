# This wrapper rescales the peaks output by macs2 to have scores between 10 and 1000

__author__ = "Chris Preusch"
__email__ = "cpreusch@ust.hk"

from snakemake.shell import shell
import os

scores_col=5

def get_min_score(fn):
    sort_cmd = "tail -n+2 {fn} | sort -k {scores_col}gr,{scores_col}gr".format(scores_col=scores_col, fn=fn)
    cmd = "{sort_cmd} | awk 'BEGIN{{{{FS=\"\\t\";OFS=\"\\t\"}}}}{{{{if (NF != 0) print $0}}}}' | tail -n 1 - | cut -f {scores_col}".format(sort_cmd=sort_cmd, scores_col=scores_col)
    ret = shell(cmd, read=True)
    return float(ret.strip())

def get_max_score(fn):
    sort_cmd = "tail -n+2 {fn} | sort -k {scores_col}g,{scores_col}g".format(scores_col=scores_col, fn=fn)
    cmd = "{sort_cmd} | awk 'BEGIN{{{{FS=\"\\t\";OFS=\"\\t\"}}}}{{{{if (NF != 0) print $0}}}}' | tail -n 1 - | cut -f {scores_col}".format(sort_cmd=sort_cmd, scores_col=scores_col)
    ret = shell(cmd, read=True)
    return float(ret.strip())

def rescale_scores_cmd(fn, new_min=10, new_max=1000):
    a = get_min_score(fn)
    b = get_max_score(fn)
    x = new_min
    y = new_max

    if a == b:
        rescale_formula = "x"
    else:
        rescale_formula = "((n-a)*(y-x)/(b-a))+x"

    cmd = "awk 'BEGIN{{OFS=\"\\t\"}}{{n=${scores_col};a={a};b={b};x={x};y={y}}}{{${scores_col}=int({rescale_formula}) ; print $0}}'".format(scores_col=scores_col, a=a, b=b, x=x, y=y, rescale_formula=rescale_formula)
    return cmd

rescale_cmd = rescale_scores_cmd(snakemake.input.bed)

shell(
    "export LC_COLLATE=C; slopBed -i {snakemake.input.bed} -g {snakemake.input.sizes} -b 0 "
    "| bedClip stdin {snakemake.input.sizes} stdout "
    "| {rescale_cmd} | sort -k1,1d -k2,2n > {snakemake.output[0]}"
    )
