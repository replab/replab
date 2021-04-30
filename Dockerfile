# Use Octave docker image
FROM       gnuoctave/octave:6.2.0

# Usage:
# docker run -it -v <your directory>:/documents/

# Install base packages
RUN apt-get update && apt-get install -y -q python3-pip pandoc unzip

RUN apt-get install build-essential
RUN apt-get install liboctave-dev


# Install Python package

RUN mkdir /workspace && mkdir /workspace/replab

RUN chmod 777 -R /workspace

COPY sphinx/requirements.txt /workspace

#The following line used to produce the error: "The command '/bin/sh -c pip3 install -r /workspace/requirements.txt' returned a non-zero code: 2"
#RUN pip3 install -r /workspace/requirements.txt

#WORKDIR /workspace/replab
VOLUME /workspace/replab

CMD ["/bin/bash"]

RUN pwd
RUN ls
RUN ls /workspace
RUN ls /workspace/replab

# Copies your code file from your action repository to the filesystem path `/` of the container
#COPY entrypoint.sh /workspace/entrypoint.sh


# Copies the whole repository
COPY . /workspace/replab/
RUN ls /workspace/replab
RUN ls /workspace/replab/external/sdpt3

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/workspace/replab/entrypoint.sh"]

