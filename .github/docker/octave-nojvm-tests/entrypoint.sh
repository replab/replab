#!/bin/bash -l

# This script runs the tests with or without covering

echo "Executing entrypoint.sh"

echo argument=$1

pwd
ls -al

git rev-parse HEAD

export ADDPATH_COMMAND="replab_init('verbose', 2);"
export COVERING=false

if [ $COVERING == true ]; then
  TEST_COMMAND="exit(~replab_runtests(1,1));";
else
  TEST_COMMAND="exit(~replab_runtests(0,1));";
fi
echo "TEST_COMMAND=$TEST_COMMAND";

# Check what octave packages we have installed
octave -q --eval "ver"

# To check that octave cannot access java
echo "The following line should give an error ---"
octave --eval "b = javaMethod('valueOf', 'java.math.BigInteger', 2)"
echo "--- end of the expected error"

# Remove any cached results files from previous build, if present
rm -f testresults.xml;

# Run tests
if octave -q --eval "$ADDPATH_COMMAND $TEST_COMMAND"; then
  # Check where we ended up and what's going on where we are
  pwd
  ls -alh
  if [ $COVERING == true ]; then
    bash <(curl -s https://codecov.io/bash);
  fi
else
  # The tests did not pass
  exit 1
fi
