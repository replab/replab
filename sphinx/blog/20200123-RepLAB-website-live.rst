The RepLAB documentation and the infrastructure behind it
=========================================================

.. post:: 26 Jan, 2020
   :tags: release
   :category: Relase
   :author: Denis Rosset
   :location: Waterloo, Ontario, Canada
   :excerpt: 2
   :image: 1
   :nocomments:

Our ambition for RepLAB is two-fold: it makes hard problems tractable by reducing the complexity using symmetries (and it is currently used in published results), and it makes an excellent tool to teach group and representation theory.

Until now, no immediately accessible tool was available for beginners in those topics. `GAP System <https://www.gap-system.org/>`_ and `SageMath <https://www.sagemath.org/>`_ are pretty incredible achievements, but they do not directly integrate with environments used in first year of undergraduate teaching. RepLAB can provide that experience, and be used to explore group representations numerically in the MATLAB/`Octave <https://www.gnu.org/software/octave/>`_ environment familiar to students.

MATLAB has comprehensive documentation, and an internal help system integrated with the REPL (the `help` command). However, MathWorks does not provide a comprehensive documentation solution such as `Sphinx <http://www.sphinx-doc.org/en/master>`_ that can integrate documentation with API reference.

We based our documentation on Sphinx through the `Sphinx MATLAB domain <https://github.com/sphinx-contrib/matlabdomain>`_, which now powers this documentation website. To complement Sphinx abilities, in particular in the support of inheritance, we resorted to code generation.

With the recent work on the RepLAB infrastructure (see `+replab.+infra.CodeBase`), we are able to parse the MATLAB source files, and process them by adding the missing elements in the API documentation comments, in particular inherited methods. Our parser is not a general MATLAB source code parser, but handles the particular subset used in RepLAB well.

.. figure:: GenerateSystem.gif
   :width: 483px
   :height: 361px
   :align: center
   :figclass: align-center

   Document generation in RepLAB through the `replab_generate` command.

The same infrastructure is used to power our replacement for the MATLAB ``help`` command. Indeed, the ``help`` command in MATLAB does not handle inheritance gracefully. We also want to give a good user experience in the open source Octave clone. We thus process the Sphinx documentation comments so that they can be displayed in the MATLAB REPL while preserving as much of the interactivity and formatting as possible.

.. figure:: ../HelpSystem.gif
   :width: 483px
   :height: 361px
   :align: center
   :figclass: align-center

   Integrated help system in RepLAB

We also switched from <MATLAB Live Scripts `https://www.mathworks.com/help/matlab/live-scripts-and-functions.html`>_ to <Jupyter `https://jupyter.org/`>_ for our tutorials, which unlocks great integration with Sphinx.

Currently, all this infrastructure runs on `Octave <https://www.gnu.org/software/octave/>`_ as well, which is important for our `continuous integration <https://travis-ci.com/replab/replab>`_ process (we still wish MathWorks had a cheap license for CI in open source projects, as Octave can be around 10x slower).
