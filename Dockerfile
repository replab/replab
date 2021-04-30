# Use Octave docker image
FROM       gnuoctave/octave:6.2.0

# Usage:
# docker run -it -v <your directory>:/documents/

# Install base packages
RUN apt-get update && apt-get install -y -q python3-pip pandoc unzip

# Install Python package

RUN mkdir /workspace && mkdir /workspace/replab

RUN chmod 777 -R /workspace

COPY sphinx/requirements.txt /workspace

RUN pip3 install -r /workspace/requirements.txt

#WORKDIR /workspace/replab
VOLUME /workspace/replab

CMD ["/bin/bash"]

# Copies the whole repository
COPY . /workspace/replab/

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/workspace/replab/entrypoint.sh"]

