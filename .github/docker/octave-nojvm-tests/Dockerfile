# Use Octave docker image
FROM       gnuoctave/octave:5.2.0

# Usage:
# docker run -it -v <your directory>:/documents/

USER root
# Install base packages
RUN apt-get update && apt-get install -y -q python3-pip pandoc unzip curl

# Remove java
RUN sudo update-alternatives --remove-all "java"
RUN sudo update-alternatives --remove-all "javac"
RUN sudo apt-get -y -q purge openjdk-\*

RUN mkdir /workspace

RUN chmod 777 -R /workspace

CMD ["/bin/bash"]

# Copies the entrypoint script
COPY entrypoint.sh /

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/entrypoint.sh"]

