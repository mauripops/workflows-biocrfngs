#!/usr/bin/env python
import os
import subprocess
import logging
import shutil
import csv
import sys
import warnings
import pandas as pd
warnings.simplefilter(action='ignore', category=DeprecationWarning)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("post_gen_project")

PROJECT_DIRECTORY = os.path.realpath(os.path.curdir)
SAMPLE_SHEET = "{{ cookiecutter.sample_sheet }}"
WORKFLOW_EXPERIMENT_ID = "{{ cookiecutter.experiment_id }}"
SPECIES = "{{ cookiecutter.ref_species.replace(' ','_') }}"

root = os.getcwd()

if os.path.exists(SAMPLE_SHEET):
    # dst = os.path.join(root, "sample_sheet.csv")
    # shutil.copyfile(sample_sheet, dst)
    sample = pd.read_csv(SAMPLE_SHEET, delimiter="\t")
    
    sample.rename(columns = {
        "Sample name": "sample",
        "read1":"fq1",
        "read2":"fq2"},
        inplace=True)
    samples_tsv = pd.Series(sample["sample"].unique())
    samples_tsv = pd.DataFrame(samples_tsv, columns=["sample"])

    sample['unit'] = sample.groupby('Sample number').cumcount()
    units_tsv = sample[["sample","unit","fq1","fq2"]]
    units_tsv.loc[:,"fq1"] = units_tsv["fq1"].apply(lambda x: os.path.join("data/reads",x))
    units_tsv.loc[:,"fq2"] = units_tsv["fq2"].apply(lambda x: os.path.join("data/reads",x))
    
    samples_tsv.to_csv(os.path.join(PROJECT_DIRECTORY, WORKFLOW_EXPERIMENT_ID,"config/samples.tsv"),
                       index=False, sep="\t")
    units_tsv.to_csv(os.path.join(PROJECT_DIRECTORY, WORKFLOW_EXPERIMENT_ID,"config/units.tsv"),
                     index=False, sep="\t")
    os.symlink('/data2/biocrfngs/reference_genomes/{}'.format(SPECIES), os.path.join(PROJECT_DIRECTORY, WORKFLOW_EXPERIMENT_ID, 'data', SPECIES))

else:
    # no sample sheet, exit
    print("ERROR: The sample_sheet given could not be found: {}".format(SAMPLE_SHEET))
    sys.exit(1)
