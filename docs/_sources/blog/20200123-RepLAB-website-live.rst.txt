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

Our ambition for RepLAB is two-fold: reduce the complexity of hard problems using their symmetries (we used it in published results), and provide as well an accessible tool to teach group and representation theory.

Until now, no immediately accessible tool was available for beginners in those topics. `GAP System <https://www.gap-system.org/>`_ and `SageMath <https://www.sagemath.org/>`_ are pretty incredible pieces of software, but they do not directly integrate with environments used in the first year of undergraduate teaching. RepLAB can provide that experience, and be used to explore group representations numerically in the `MATLAB <https://www.mathworks.com/products/matlab.html>`_/`Octave <https://www.gnu.org/software/octave/>`_ environment familiar to students.

It was thus paramount for us to document the use our code well.

MATLAB has comprehensive documentation, and an internal help system integrated with the REPL (the `help` command). Still, MathWorks does not provide a comprehensive documentation solution such as `Sphinx <http://www.sphinx-doc.org/en/master>`_ that can integrate documentation with API reference and publish a website.

We based our documentation on Sphinx using the `Sphinx MATLAB domain <https://github.com/sphinx-contrib/matlabdomain>`_, and it now powers this documentation website. To complement Sphinx abilities, in particular in the support of inheritance, we resorted to code generation.

With the recent work on the RepLAB infrastructure (see `+replab.+infra.CodeBase`), we are able to parse the RepLAB source files, and process them by adding the missing elements in the API documentation comments, in particular inherited methods. Our parser is not a general MATLAB source code parser, but handles the particular subset used in RepLAB well.

.. figure:: GenerateSystem.gif
   :align: center
   :figclass: align-center

   Document generation in RepLAB through the `root.replab_generate` command.

The same infrastructure is used to power our replacement for the MATLAB ``help`` command (for example, the ``help`` command in MATLAB does not handle inheritance gracefully). We also want to give a good user experience in the open source Octave clone (and the `help` command in Octave does not handle much object-oriented programming). We thus take charge of processing the Sphinx documentation comments so that they can be displayed in the MATLAB REPL while preserving as much of the interactivity and formatting as possible.

.. figure:: ../HelpSystem.gif
   :width: 483px
   :height: 361px
   :align: center
   :figclass: align-center

   Our integrated help system in RepLAB

In the website construction, we also switched from `MATLAB Live Scripts <https://www.mathworks.com/help/matlab/live-scripts-and-functions.html>`_ to `Jupyter <https://jupyter.org/>`_ for our tutorials, as they integrate beautifully with Sphinx.

Currently, all this infrastructure runs on `Octave <https://www.gnu.org/software/octave/>`_ as well, which is important for our `continuous integration <https://travis-ci.com/replab/replab>`_ process. We still wish MathWorks had a cheap license for CI in open source projects, as Octave can be around 10x slower!

We are now focusing on improving the RepLAB core algorithms, in particular by performing proper analysis of the numerical errors, while enriching the documentation with guided tutorials. Stay tuned!
