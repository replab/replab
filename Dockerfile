# Use Ubuntu LTS release
FROM       ubuntu:18.04

# Usage:
# docker run -it -v <your directory>:/documents/

# Install base packages
RUN apt-get update && apt-get install -y -q python3-pip pandoc octave octave-optim

# Install Python package

RUN pip3 install sphinx sphinxcontrib-matlabdomain texext guzzle_sphinx_theme sphinx-collapse-admonitions sphinxcontrib-fulltoc ablog jupyter_client octave-kernel nbsphinx

RUN mkdir /workspace && mkdir /workspace/replab

WORKDIR /workspace/replab
VOLUME /workspace/replab

CMD ["/bin/bash"]
