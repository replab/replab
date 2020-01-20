#!/usr/bin/env python

import shutil
import subprocess
import sys

# Creates static website for the API documentation, and places it in
# the right location for jekyll.

# clear output directory
shutil.rmtree('_build', True);

# runs Sphinx
status = subprocess.call(['make', 'html']);
if status != 0:
  sys.exit(status)

# copy folder to jekyll
status = shutil.copytree('_build', '../docs/sphinx')
if status != 0:
  sys.exit(status)

