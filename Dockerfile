# Use Ubuntu LTS release
FROM       ubuntu:18.04

# Usage:
# docker run -it -v <your directory>:/documents/

# Install base packages
RUN apt-get update && apt-get install -y -q python3-pip pandoc octave octave-optim liboctave-dev openjdk-11-jre

# Install Python package

RUN mkdir /workspace && mkdir /workspace/replab

RUN chmod 777 -R /workspace

COPY sphinx/requirements.txt /workspace

#RUN pip3 install -r /workspace/requirements.txt

WORKDIR /workspace/replab
VOLUME /workspace/replab

CMD ["/bin/bash"]

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/workspace/replab/entrypoint.sh"]

