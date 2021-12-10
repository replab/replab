.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials/index
   howto/index
   topic/index
   publications/index
   api/index
   development/index
   blog/index

Welcome to RepLAB!
==================

Current version: |version| (`GitHub <https://github.com/replab/replab>`_ / `latest release ZIP`_ / :doc:`installation instructions <tutorials/installation>`).

**RepLAB** provides tools to study representations of finite groups and decompose them numerically.
It is compatible with both `MATLAB <https://www.mathworks.com/products/matlab.html>`_ and `Octave <https://www.gnu.org/software/octave/>`_.

Key features
------------

.. csv-table::
   :header: "Provable decompositions", " ", "Modular construction", " ", "Optimization toolbox integration", " ", "Interactive documentation", " ", "Open source"
   :header-rows: 1
   :widths: 20, 1, 20, 1, 20, 1, 20, 1, 20
   
   " ", " ", " ", " ", " ", " ", " ", " ", " "
   " ", " ", " ", " ", " ", " ", " ", " ", " "
   "RepLAB's numerical decomposition into irreps can be certified", " ", "Construct new groups and representations by simple combination", " ", "Solves convex optimization problems with symmetries efficiently", " ", "Extensive documentation accessible by a click", " ", "Full open source software compatible with the GNU/Octave interpreter"

RepLAB in action
----------------

.. figure:: CompactGroups.gif
   :width: 483px
   :align: center
   :figclass: align-center responsive

   Decomposition of the :math:`U \otimes U` representation of the unitary group of dimension 2.


.. figure:: HelpSystem.gif
   :width: 483px
   :align: center
   :figclass: align-center responsive

   Integrated help system in RepLAB


How to start using RepLAB?
--------------------------

Have a look at the :doc:`tutorials <tutorials/index>`!
They have everything to get you started, from installation instructions to hands-on examples.

The documentation of **RepLAB** is organized along 4 directions, following this `approach <https://www.divio.com/blog/documentation/>`_:

-  :doc:`Tutorials <tutorials/index>`: are short hands-on presentations that give you a taste of the goodness of **RepLAB**
-  :doc:`How-to guides <howto/index>`: are concise recipes showing how to achieve a specific goal
-  :doc:`Topic guides <topic/index>`: are understanding-oriented presentations explaining the big picture and key notions on which this software is built
-  :doc:`Technical reference <api/index>`: contains a complete and accurate description of each object of the library

Work in progress
----------------

**RepLAB** is a work-in-progress. In particular:

- We are still working on the estimation/control of numerical errors (a few criteria are currently hard-coded).
  That said, **RepLAB** works fine for representations of medium size ($d$ equal to a few thousands).
- While RepLAB implements several computational group theory algorithms, their performance is not competitive with dedicated computational group theory software such as GAP or Magma.
- The documentation is still a work in progress, and tutorials/how to guides/... will be added as the project progresses.

Why RepLAB?
-----------

Because no open source library decomposes arbitrary permutation/monomial representations into irreducible representations over the reals. **RepLAB** implements numerical methods that perform this decomposition up to machine precision.

That said, other libraries working on the same problem space include:

-  The `GAP System 4 package RepnDecomp <https://gap-packages.github.io/RepnDecomp>`_ by
   Kaashif Hymabaccus.

-  The `GAP System 3 package AREP <https://www.gap-system.org/Gap3/Packages3/arep.html>`_ by
   Sebastian Egner and Markus Püschel.

-  `NCSOStools <http://ncsostools.fis.unm.si/documentation/awbd>`_ includes an implementation of the Murota-Kanno-Kojima-Kojima-Maehara algorithm to decompose matrix \*-algebras.

Additional features unique to **RepLAB**:

-  **RepLAB** follows a category-based approach. This enables preserving the structure of objects under operations, naturally allowing composition of finite and continuous groups/representations together for instance.
-  **RepLAB** integrates with the convex optimization framework `YALMIP <https://github.com/yalmip/YALMIP>`_. This allows to formulate optimizations subject to symmetries easily.


Contributors
------------

**RepLAB** and the group theory/linear algebra libraries it depends on were initiated by `Denis Rosset <https://github.com/denisrosset>`_ and `Jean-Daniel Bancal <https://github.com/jdbancal>`_. The project now has more `contributors <https://github.com/replab/replab/graphs/contributors>`_.

**RepLAB** references in the ``/external`` directory the following libraries: the `MOxUnit <https://github.com/MOxUnit/MOxUnit>`_ test framework by Nikolaas N. Oosterhof, the `YALMIP <https://github.com/yalmip/YALMIP>`_ toolbox for optimization modeling by Johan Löfberg, the `SDPT3 <https://github.com/sqlp/sdpt3>`_ solver, and the `VPI <https://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic>`_ big integer library by John D'Errico.

Feedback and suggestions are always welcome. We ask participants to follow the guidelines of the `Typelevel Code of Conduct <https://typelevel.org/conduct.html>`_.

License
-------

**RepLAB** is (C) 2018-2021 Denis Rosset, Jean-Daniel Bancal and other collaborators, and licensed under the `Mozilla Public License 2.0 <https://github.com/replab/replab/LICENSE>`_.
