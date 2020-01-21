Sphinx tips and tricks
======================

We write here miscellaneous tips and tricks.

Sphinx and Jupyter required packages
------------------------------------

I tried with Python 3 in a dedicated virtualenv in Ubuntu 18.04 as follows:

Run somewhere

::

    python3 -m venv sphinx

Then enter the environment

::

    source sphinx/bin/activate

Install the following python packages:

- sphinx
- sphinxcontrib-matlabdomain
- texext
- guzzle_sphinx_theme
- sphinx-collapse-admonitions
- sphinxcontrib-fulltoc

This can be done with the command

::

   pip install sphinx sphinxcontrib-matlabdomain texext guzzle_sphinx_theme sphinx-collapse-admonitions sphinxcontrib-fulltoc

(add the option "--user" if not installing inside a virtual environment)

For the Jupyter notebooks:

Follow `<https://nbsphinx.readthedocs.io/en/latest/installation.html>`_

In particular, install the Python package nbsphinx, and the Ubuntu package pandoc


Later on?

Proof environment? `<https://framagit.org/spalax/sphinxcontrib-proof/>`_

Resources: `<https://github.com/yoloseem/awesome-sphinxdoc>`_

References to root-level scripts
--------------------------------

Scripts such as `~root.replab_init`, `~root.replab_generate` live in the fake ``root`` package, so these references look actually like `root.replab_init` and `root.replab_generate`.

The ``root`` prefix is stripped when displayed on the console. In the Sphinx documentation, one should use the ``~`` prefix.

To convert Matlab Live Scripts to Jupyter notebooks
---------------------------------------------------

`Jupytext <https://github.com/mwouts/jupytext>`_ works fantastically well.

Generation of Jupyter notebooks
-------------------------------

Running the Matlab kernel from within a Matlab session does not work. We sidestep this problem by using the Octave kernel in all Jupyter notebooks.
