Docker Octave image
===================

Due to the fragility of the combination of the Sphinx/Matlab/Octave ecosystems, we created a Docker image able to run the tests and compile the documentation.

The latest image can be pulled using the command:

::
   docker pull replab/replab:latest

and the image can be run using the command:

::
   docker run --user $(id -u):$(id -g) -e HOME=/workspace -v /PATH/TO/LOCAL/REPLAB/REPO:/workspace/replab -it replab/replab:latest


The image is created using the ``root`` user, albeit giving full permissions to others inside the container; thus additional files can be created by the appropriate non-root user when running it.
