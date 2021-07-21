Named discrete groups
=====================

.. module:: +replab

RepLAB provides constructions of standard discrete groups.

General (signed) permutation groups are the two groups below.

* `.SymmetricGroup`: Describes the symmetric group on a given domain.
* `.S`: Shorthand for the construction of a symmetric group.
* `.SignedSymmetricGroup`: Describes the signed symmetric group.

Then we distinguish particular groups.

* `.A`: Describes the alternating group on a given domain.
* `.C`: Describes the cyclic group on a given domain.
* `.D`: Describes the dihedral group on a given domain.
* `.KleinFourGroup`: Describes the Klein four-group
* `.QuaternionGroup`: Describes the quaternion group of order 8.

SymmetricGroup
++++++++++++++

.. autoclass:: SymmetricGroup

S
+

.. autofunction:: S

SignedSymmetricGroup
++++++++++++++++++++

.. autoclass:: SignedSymmetricGroup

A
+

.. autofunction:: A

C
+

.. autofunction:: C

D
+

.. autofunction:: D


KleinFourGroup
++++++++++++++

.. autofunction:: KleinFourGroup

QuaternionGroup
+++++++++++++++

.. autofunction:: QuaternionGroup
