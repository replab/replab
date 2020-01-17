Current version: **{{site.replabVersion}}**.

|Join the chat at https://gitter.im/denisrosset/replab| |Travis CI|
|codecov|

**RepLAB** provides tools to study representations of finite groups and
decompose them numerically. It is compatible with both MATLAB and
Octave.

How to start using RepLAB?
--------------------------

Have a look at the `**tutorials** <docs/tutorials/tutorials.html>`__!
They have everything to get you started, from installation instructions
to hands-on examples.

The documentation of **RepLAB** is organized along 4 directions,
following a recent `discussion on software
documentation <https://www.divio.com/blog/documentation/>`__:

-  `Tutorial <docs/tutorials/tutorials.html>`__: are short hands-on
   presentations that give you a taste of the goodness of **RepLAB**
-  `How-to guides <docs/howto/howto.html>`__: are concise recipes that
   show you how to achieve a specific goal
-  `Topic guides <docs/topic/guides.html>`__: are understanding-oriented
   presentations that explain the big picture and the key notions on
   which this software is built
-  `Technical reference <docs/reference/reference.html>`__: contains a
   complete and accurate description of each object of the library

Work in progress
----------------

**RepLAB** is a work-in-progress. In particular:

-  RepLAB only works in double floating-point precision.
-  Most representations are handled by decomposing preimages in words
   over generators, which limits the order of groups that can be
   handled.
-  While RepLAB has a basic implementation of the BSGS construction, it
   does not offer much to work with permutation groups.

Why RepLAB?
-----------

Because no open source library exists to decomposes arbitrary
permutation/monomial representations into irreducible representations
over the reals. RepLAB implements numerical methods that perform this
decomposition up to machine precision.

That said, other libraries working on the same problem space include:

-  The `GAP System 3 package
   AREP <https://www.gap-system.org/Gap3/Packages3/arep.html>`__ by
   Sebastian Egner and Markus Püschel.
-  `NCSOStools <http://ncsostools.fis.unm.si/documentation/awbd>`__
   includes an implementation of the Murota-Kanno-Kojima-Kojima-Maehara
   algorithm to decompose matrix \*-algebras.

Documentation and Support
-------------------------

-  Chat it up on `Gitter <https://gitter.im/denisrosset/replab>`__.
-  Check the `tutorial <docs/installation.html>`__.

Contributors
------------

RepLAB and the group theory/linear algebra libraries it depends on were
initiated by `Denis Rosset <https://github.com/denisrosset>`__ and
`Jean-Daniel Bancal <https://github.com/jdbancal>`__. The project has
now more
`contributors <https://github.com/replab/replab/graphs/contributors>`__.

RepLAB references in the ``/external`` directory the following
libraries:

-  As a Git submodule, the
   `MOxUnit <https://github.com/MOxUnit/MOxUnit>`__ test framework by
   Nikolaas N. oosterhof.

-  A copy of the
   `VPI <https://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic>`__
   big integer library by John D'Errico.

Feedback and suggestions are always welcome. We ask participants to
follow the guidelines of the `Typelevel Code of
Conduct <https://typelevel.org/conduct.html>`__.

License
-------

RepLAB is (C) 2018-2019 Denis Rosset, Jean-Daniel Bancal and other
collaborators, and licensed under the `Mozilla Public License
2.0 <https://github.com/replab/replab/LICENSE>`__.

.. |Join the chat at https://gitter.im/denisrosset/replab| image:: https://badges.gitter.im/Join%20Chat.svg
   :target: https://gitter.im/denisrosset/replab
.. |Travis CI| image:: https://travis-ci.com/replab/replab.svg?branch=master
   :target: https://travis-ci.com/replab/replab
.. |codecov| image:: https://codecov.io/gh/replab/replab/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/replab/replab
