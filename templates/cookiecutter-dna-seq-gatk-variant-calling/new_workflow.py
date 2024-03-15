#!/usr/bin/env python

import os
from cookiecutter.main import cookiecutter

# Create project from the template in this folder
cookiecutter(os.path.dirname(os.path.realpath(__file__)),
    overwrite_if_exists=True,
    no_input=True
    )
