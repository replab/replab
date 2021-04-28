#!/bin/bash -l

echo "Executing entrypoint.sh"

export ADDPATH_COMMAND="replab_init('verbose', 2);"
export COVERING=false

pwd
ls

if [ $COVERING == true ]; then
  TEST_COMMAND="exit(~replab_runtests(1,1));";
else
  TEST_COMMAND="exit(~replab_runtests(0,1));";
fi
echo "TEST_COMMAND=$TEST_COMMAND";

# Double-check we are still in the right directory
pwd

# Check what octave packages we have installed
octave -q --eval "ver"
# Check that octave can access java
octave --eval "b = javaMethod('valueOf', 'java.math.BigInteger', 2)"

# Remove any cached results files from previous build, if present
rm -f testresults.xml;

# Run tests
octave -q --eval "$ADDPATH_COMMAND $TEST_COMMAND";

# Check where we ended up and what's going on where we are
pwd
ls -alh
if [ $COVERING == true ]; then
  bash <(curl -s https://codecov.io/bash);
fi

