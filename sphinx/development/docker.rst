Docker Octave image
===================

Due to the fragility of the combination of the Sphinx/Matlab/Octave ecosystems, we created a `Docker <https://docs.docker.com/install/>`_ image able to run the tests and compile the documentation.

The latest image can be pulled using the command:

::

   docker pull replab/replab:latest

and the image can be run using the command:

::

   docker run --user $(id -u):$(id -g) -e HOME=/workspace -v /PATH/TO/LOCAL/REPLAB/REPO:/workspace/replab -it replab/replab:latest


The image is created using the ``root`` user, albeit giving full permissions to others inside the container; thus additional files can be created by the appropriate non-root user when running it.

The above command mounts the replab folder of the host machine into the container. Any change to this folder inside the container will thus be immediately reflected in the host filesystem.


Update dependencies versions in the Docker file
-----------------------------------------------

The underlying Linux distribution is Ubuntu 18.04 LTS, so we can expect some stability from the packages managed by ``apt``. However, we pin the Python packages versions using `pip-tools <https://pythonspeed.com/articles/pipenv-docker/>`_.

In the ``sphinx`` subdirectory, the file ``requirements.in`` lists the packages we require, while the file ``requirements.txt`` lists the python packages we require. Assuming that the latest dependencies form consistent python environment able to compile the documentation, the ``requirements.txt`` file can be updated with the command:

::

   pip-compile requirements.in > requirements.txt

run from the ``sphinx/`` directory.
