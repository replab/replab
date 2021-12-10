Installation
============

In this chapter we cover the basics of starting to use **RepLAB**.

Downloading the library
-----------------------

The representation theory code of the library is self-contained.
Extended features such as *unit tests*, *code coverage* and *convex optimization*
make use of external code. Here are two ways of installing the library with the
desired set of features to get started. Choose the one which suits you best.

Option 1: Download the latest RepLAB release, and use our easy install script.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download the `latest release ZIP`_, see the `GitHub release page <https://www.github.com/replab/replab/releases>`_ for all releases.
   This ZIP file includes the core code, but not the external dependencies.

2. Launch MATLAB/Octave, run the ``replab_init autoinstall`` script in the root folder
   which will take care of downloading and installing dependencies, which are:

-  to run tests: `MOxUnit <https://github.com/MOxUnit/MOxUnit>`__
-  to define convex optimization (SDP) problems and run corresponding
   tests: `YALMIP <https://github.com/yalmip/YALMIP>`__
-  to solve SDP problems and run corresponding tests:
   `SDPT3 <https://github.com/sqlp/sdpt3>`__

3. Go to the "Initializing the library" section below and follow the usage instructions there.

..
   .. figure:: EasyInstall.gif
   :width: 483px
   :height: 361px
   :align: center
   :figclass: align-center

   From the release ZIP file to a successful RepLAB installation.


Option 2: Clone the library
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: For advanced users only
   :class: collapsed

   Clone the library from GitHub using the following command:

   ::

       git clone --recursive https://www.github.com/replab/replab

   which will download the latest ``master`` version, and update the Git
   submodules automatically.

   This creates a folder **RepLAB** with all the necessary code, including the
   `VPI <https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic>`__
   library for large integers, the testing suite
   `MOxUnit <https://github.com/MOxUnit/MOxUnit>`__, and the tools needed
   for semidefinite programming.

   To ignore changes that will happen when the SDPT3 submodule is compiled,
   (and thus untracked files are created), run the following command:

   ::

      git config submodule.SDPT3.ignore untracked

   or add a ``ignore = untracked`` line to the ``[submodule "SDPT3"]`` section
   in ``.git/config``.

   If you want to run the release process, check out the ``gh-pages`` branch in the
   ``docs/`` folder by running:

   ::

      git worktree add docs gh-pages

   in the root RepLAB folder.


Initializing the library
------------------------

To use the library, the **RepLAB** folder must be added in Matlab or Octave.
Additional paths are also necessary to enable specific functionalities, as
mentioned above, and a few variables mush be initialized. This can be done with

::

    replab_init

which should lead to the following output

::

    >> replab_init
    Adding RepLAB to the path
    Initializing dependency vpi
    Initializing dependency YALMIP
    Initializing dependency sdpt3
    Adding embedded SDPT3 solver to the path
    Initializing dependency MOxUnit
    >>

This command checks in particular whether an instance of YALMIP is
`available <https://yalmip.github.io/download/>`__ and
`configured <https://yalmip.github.io/tutorial/installation/>`__ on your
system. If this is not the case, the embedded version of yalmip is used.
**RepLAB** uses the `YALMIP <https://yalmip.github.io>`__ interface to
solve convex optimization problems. The ``replab_init`` command also
ensures that an `SDP solver <https://yalmip.github.io/allsolvers/>`__
such as `SeDuMi <https://github.com/SQLP/SeDuMi>`__ is properly set up.
If this is not the case, it activates the embedded SDPT3 solver. The
proper installation of a YALMIP instance can be checked with the command
``yalmiptest``.

The command ``replab_init`` should always be used before running any
**RepLAB** command. This command only takes some time to run the first
time it is called in a MATLAB/Octave session.

Customizing the initialization script
-------------------------------------

The script ``replab_config.m`` present in the root folder can be customized:
for example, unneeded dependencies can be commented out, and various RepLAB
options can be configured.

This script is run once during the ``replab_init`` initialization; subsequent calls
to ``replab_init`` will not run it again, except if a ``clear all`` command has been
run in between.

Testing
-------

The proper installation of **RepLAB** can then be checked by running the
test commands:

::

    replab_runtests

This checks the proper working of the whole package (requires the
companion packages for test and convex optimization).
