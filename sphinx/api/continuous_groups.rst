Continuous groups
=================

.. module:: +replab

**RepLAB** provides constructions of standard continuous groups.

First, we consider torus-like groups.

* `.TorusGroup`: Describes a compact, connected, abelian Lie group (i.e. a torus group).
* `.T`: Describes the standard torus group of given rank.

Then, we have the classical compact groups.

* `.ClassicalCompactGroup`: Class implementing all the groups below.
* `.O`: Constructs the group of (real) orthogonal matrices of a given size.
* `.SO`: Constructs the group of (real) special orthogonal matrices of a given size.
* `.U`: Constructs the group of (complex) unitary matrices of a given size.
* `.SU`: Constructs the group of (complex) special unitary matrices of a given size.
* `.Sp`: Constructs the compact symplectic group of a given size, which encodes using complex numbers the generalization of the unitary group to quaternions.

TorusGroup
++++++++++

.. autoclass:: TorusGroup

T
+

.. autofunction:: T

ClassicalCompactGroup
+++++++++++++++++++++

.. autoclass:: ClassicalCompactGroup

O
+

.. autofunction:: O

SO
++

.. autofunction:: SO

U
+

.. autofunction:: U

SU
++

.. autofunction:: SU

Sp
++

.. autofunction:: Sp
