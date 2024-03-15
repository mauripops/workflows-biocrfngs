#!/usr/bin/env python
import os
import subprocess
import logging
import shutil
import csv
import sys
import pandas as pd

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("pre_gen_project")
