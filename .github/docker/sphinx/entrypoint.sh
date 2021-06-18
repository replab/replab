#!/bin/bash -l

# This script builds the sphinx documentation

echo "Executing entrypoint_sphinx.sh"

echo argument=$1

pwd
ls -al

git rev-parse HEAD

# Prepare doc folder and install python requirements
ls -al docs
pip3 install -r sphinx/requirements.txt

# prepare commands
export ADDPATH_COMMAND="replab_init('verbose', 2);"
GENERATE_COMMAND="exit(~replab_generate('sphinx'));";
echo "GENERATE_COMMAND=$GENERATE_COMMAND";

# Check what octave packages we have installed
octave -q --eval "ver"

# Check that octave can access java
octave --eval "b = javaMethod('valueOf', 'java.math.BigInteger', 2)"

# Run commands
if octave -q --eval "$ADDPATH_COMMAND $GENERATE_COMMAND"; then
  # Check where we ended up and what's going on where we are
  pwd
  ls -alh
else
  # The commands failed
  exit 1
fi
