#!/usr/bin/env python3
import sys
import re
import argparse
import subprocess
import os

from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(add_help=True)
# parser.add_argument(
#     "--help", help="Display help message.", action="store_true")
parser.add_argument(
    "script", metavar="SCRIPT",
    help="the script to pass to sbatch")

# A subset of SLURM-specific arguments
slurm_parser = parser.add_argument_group("slurm-specific arguments")
slurm_parser.add_argument(
    "-a", "--array", help="job array index values")
slurm_parser.add_argument(
    "-A", "--account", help="charge job to specified account")
slurm_parser.add_argument(
    "--begin", help="defer job until HH:MM MM/DD/YY")
slurm_parser.add_argument(
    "-c", "--cpus-per-task", help="number of cpus required per task", default=4)
slurm_parser.add_argument(
    "-d", "--dependency",
    help="defer job until condition on jobid is satisfied")
slurm_parser.add_argument(
    "-D", "--workdir", help="set working directory for batch script")
slurm_parser.add_argument(
    "-e", "--error", help="file for batch script's standard error")
slurm_parser.add_argument(
    "-J", "--job-name", help="name of job")
slurm_parser.add_argument(
    "--mail-type", help="notify on state change: BEGIN, END, FAIL or ALL")
slurm_parser.add_argument(
    "--mail-user", help="who to send email notification for job state changes")
slurm_parser.add_argument(
    "-n", "--ntasks",
    type=int,
    help="number of tasks to run (default: %(default)s)",
    default=1)
slurm_parser.add_argument(
    "-N", "--nodes",
    help="number of nodes on which to run (N = min[-max]) (default: %(default)s)",
    default=1)
slurm_parser.add_argument(
    "-o", "--output", help="file for batch script's standard output")
slurm_parser.add_argument(
    "-p", "--partition", help="partition requested")
slurm_parser.add_argument(
    "-Q", "--quiet", help="quiet mode (suppress informational messages)")
slurm_parser.add_argument(
    "-t", "--time", help="time limit")
slurm_parser.add_argument(
    "-C", "--constraint", help="specify a list of constraints")

args = parser.parse_args()

# get the last argument which should be the shell script to submit
jobscript = args.script

# get the properties from the job script
job_properties = read_job_properties(jobscript)

arg_dict = dict(args.__dict__)

job_name = "{0}.snakejob".format(job_properties["rule"])
arg_dict["job-name"] = job_name

if not os.path.exists(os.path.join(os.getcwd(), "logs")):
    print("No log directory detected, creating it.", file=sys.stderr)
    os.makedirs(os.path.join(os.getcwd(), "logs"))

arg_dict["output"] = "logs/{0}.%j.out".format(job_name)

# Process resources
if "resources" in job_properties:
    resources = job_properties["resources"]
    if arg_dict["time"] is None:
        if "runtime" in resources:
            arg_dict["time"] = resources["runtime"]
        elif "walltime" in resources:
            arg_dict["time"] = resources["runtime"]
    if "mem_per_cpu" in resources:
        arg_dict["mem-per-cpu"] = resources["mem_per_cpu"]	
        default_mem = arg_dict.pop("mem")
    elif "mem_mb" in resources:
        if resources["mem_mb"] >= 380000:
            resources["mem_mb"] = 20000
            arg_dict["mem"] = 20000
        else:
            arg_dict["mem"] = resources["mem_mb"]
    if "nodes" in resources:
        arg_dict["nodes"] = resources["nodes"]
    if "ntasks" in resources:
        arg_dict["ntasks"] = resources["ntasks"]
    gres_list = []
    if "io" in resources:
        gres_list.append("nfs:{}".format(str(resources["io"])))
    if "download" in resources:
        gres_list.append("download:{0}".format(str(resources["download"])))
    if len(gres_list) > 0:
        arg_dict["gres"] = ",".join(gres_list)

# Threads
if "threads" in job_properties:
    arg_dict["cpus-per-task"] = job_properties["threads"]
else:
    arg_dict["cpus-per-task"] = 4

# Get cluster config values
if "cluster" in job_properties:
    cluster = job_properties["cluster"]
    if arg_dict["partition"] is None:
        if "partition" in cluster:
            arg_dict["partition"] = cluster["partition"]

opt_keys = ["array", "account", "begin", "cpus-per-task",
            "depedency", "workdir", "error", "job-name", "mail-type",
            "mail-user", "ntasks", "nodes", "output", "partition",
            "quiet", "time", "constraint", "mem", "mem-per-cpu","gres"]

opt_strs = []
for k, v in arg_dict.items():
    if k not in opt_keys:
        continue
    if v is not None:
        opt_strs.append("--{}=\"{}\"".format(k, v))

cmd = "sbatch --begin=now+3 {opts} {script}".format(opts=" ".join(opt_strs), script=args.script)

# print the sbatch command to stderr
print(cmd, file=sys.stderr)

# Used for testing
# print(cmd)
# sys.exit(0)

# Try to run the command
try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise
