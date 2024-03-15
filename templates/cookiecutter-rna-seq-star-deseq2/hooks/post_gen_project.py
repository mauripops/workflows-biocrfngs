#!/usr/bin/env python
import os
import subprocess
import logging
import shutil
import csv
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("post_gen_project")

PROJECT_DIRECTORY = os.path.realpath(os.path.curdir)
SEQUENCING_FOLDER = "{{ cookiecutter.sequencing_folder }}"
WORKFLOW_USER_ID = int("{{ cookiecutter.user_id }}")
WORKFLOW_EXPERIMENT_ID = "{{ cookiecutter.experiment_id }}"
SAMPLE_SHEET = "{{ cookiecutter.sample_sheet }}"
ORGANISM = "{{ cookiecutter.ref_organism.replace(' ','_') }}"

root = os.getcwd()

# write data to filename as tsv
# data should be a list of dictionary items
def write_tsv(filename, data):
    with open(os.path.join(PROJECT_DIRECTORY, WORKFLOW_EXPERIMENT_ID, filename), 'wt') as of:
        tsv_writer = csv.writer(of, delimiter='\t')
        keys = data[0].keys()
        tsv_writer.writerow(keys)
        for d in data:
            tsv_writer.writerow([d[key] for key in keys])

if os.path.exists(SAMPLE_SHEET):
    # dst = os.path.join(root, "sample_sheet.csv")
    # shutil.copyfile(sample_sheet, dst)
    samples = pd.read_table(SAMPLE_SHEET, dtype=str).set_index(["Sample number", "Sample name", "Expt grouping"], drop=False)
    users = dict()

    zip_samples =  zip(samples["Expt grouping"].tolist(),
                       samples["Lane"].tolist(),
                       samples["Barcode"].tolist(),
                       [s.strip().replace(" ", "_") for s in samples["Sample name"].tolist()]
                       )

    split_samples = [[*sample.split("_"), lane, barcode, sample_name] for sample,lane,barcode,sample_name in zip_samples]
    files_to_download = []
    sample_tsv = []
    unit_tsv = []

    for user_id, expt_reptype, rep, lane_str, barcode, sample_name in split_samples:
        user_id = int(user_id)
        experiment_id = expt_reptype[0]

        # skip this sample if it doesn't belong to this user or experiment
        if user_id != WORKFLOW_USER_ID or experiment_id != WORKFLOW_EXPERIMENT_ID:
            continue

        condition = int(expt_reptype[1])
        rep = int(rep)

        if "," in lane_str:
            lanes = ["L0{}".format(l) for l in lane_str.split(",")]
        else:
            lanes = ["L0{}".format(lane_str)]

        barcode = barcode

        for lane in lanes:
            fq1 = "{SEQUENCING_FOLDER}_{lane}_{barcode}_1.fq.gz".format(
                        SEQUENCING_FOLDER=SEQUENCING_FOLDER,
                        lane=lane,
                        barcode=barcode)
            fq2 = "{SEQUENCING_FOLDER}_{lane}_{barcode}_2.fq.gz".format(
                        SEQUENCING_FOLDER=SEQUENCING_FOLDER,
                        lane=lane,
                        barcode=barcode)

            files_to_download.append(fq1)
            files_to_download.append(fq2)

            unit_tsv.append(dict(
                    sample=sample_name,
                    unit=lane,
                    fq1=fq1,
                    fq2=fq2
                ))

        sample_tsv.append(dict(
                sample=sample_name,
                condition="Treatment_{}".format(condition) if condition > 0 else "Control"
            ))

    if len(files_to_download) <= 0:
        # no files in analysis, something went wrong
        print("ERROR: Did not find any fastq files matching the user_id ({}) and experiment_id ({}) that you entered. Check that your sample sheet TSV file has this user and experiment, or enter a valid user_id and experiment_id when running cookiecutter.".format(user_id, experiment_id))
        sys.exit(1)
    else:
        write_tsv('samples.tsv', sample_tsv)
        write_tsv('units.tsv', unit_tsv)
        os.symlink('/data2/biocrfngs/reference_genomes/{}'.format(ORGANISM), os.path.join(PROJECT_DIRECTORY, WORKFLOW_EXPERIMENT_ID, 'resources', ORGANISM))

else:
    # no sample sheet, exit
    print("ERROR: The sample_sheet given could not be found: {}".format(SAMPLE_SHEET))
    sys.exit(1)
