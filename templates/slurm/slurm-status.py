#!/usr/bin/env python3
import re
import subprocess as sp
import shlex
import sys
import time
import logging
logger = logging.getLogger("__name__")

# constant strings that are recongnized by snakemake that give the job status
JOB_SUCCESS = 'success'
JOB_FAILED = 'failed'
JOB_RUNNING = 'running'

# try getting the status 20 times before failing
STATUS_ATTEMPTS = 20

jobid = sys.argv[1]

for i in range(STATUS_ATTEMPTS):
    try:
        sacct_res = sp.check_output(shlex.split("sacct -P -b -j {} -n".format(jobid)))
        res = {x.split("|")[0]: x.split("|")[1] for x in sacct_res.decode().split("\n")}
        break
    except sp.CalledProcessError as e:
        logger.error("sacct process error")
        logger.error(e)
    except IndexError as e:
        pass
    # Try getting job with scontrol instead in case sacct is misconfigured
    try:
        sctrl_res = sp.check_output(shlex.split("scontrol -o show job {}".format(jobid)))
        m = re.search("JobState=(\w+)", sctrl_res.decode())
        res = {jobid: m.group(1)}
        break
    except sp.CalledProcessError as e:
        logger.error("scontrol process error")
        logger.error(e)
        if i >= STATUS_ATTEMPTS - 1:
            print(JOB_FAILED)
            exit(0)
        else:
            time.sleep(1)


fail_strs = ['BOOT_FAIL', 'CANCELLED', 'FAILED', 'DEADLINE',
             'NODE_FAIL', 'PREEMPTED', 'TIMEOUT']

# Unclear whether SUSPENDED should be treated as running or failed
# fail_strs = ['BOOT_FAIL', 'CANCELLED', 'FAILED', 'DEADLINE',
#              'NODE_FAIL', 'PREEMPTED', 'TIMEOUT', 'SUSPENDED']

success_strs = ['COMPLETED']

status = res[jobid]

if status in fail_strs:
    print(JOB_FAILED)
elif status in success_strs:
    print(JOB_SUCCESS)
else:
    print(JOB_RUNNING)
