Group subsets
=============

.. module:: +replab

RepLAB has `.FiniteSet` as a base for all sets of elements coming originally from a finite group. Note that a `.FiniteGroup` itself is also a `.FiniteSet`.

* `.FiniteSet`: Base class for all subsets of a finite group.

Such subsets have a decomposition which speeds up many algorithms in RepLAB.

* `.SetProduct`: Decomposition of a `.FiniteSet` as a cartesian product of sets.

A finite group is partitioned into conjugacy classes.

* `.ConjugacyClasses`: Set of all conjugacy classes of a finite group.
* `.ConjugacyClass`: Conjugacy class of a finite group.

A finite group can also be partitioned into cosets by one (or two!) of its subgroups.

* `.Cosets`: Base class for sets of cosets of a finite group.
* `.Coset`: Base class for all cosets.

* `.CosetTable`: Describes the left action of a group on its cosets by a subgroup.

* `.LeftCosets`, `.RightCosets`, `.NormalCosets`: Sets of left/right/normal cosets, obtained using a subgroup of a finite group

* `.LeftCoset`, `.RightCoset`, `.NormalCoset`: Describes an element of the above

* `.DoubleCosets`: Sets of double cosets, obtained using two subgroups of a finite group
* `.DoubleCoset`: Describes an element of the above

FiniteSet
+++++++++

.. autoclass:: FiniteSet

SetProduct
++++++++++

.. autoclass:: SetProduct

ConjugacyClasses
++++++++++++++++

.. autoclass:: ConjugacyClasses

ConjugacyClass
++++++++++++++

.. autoclass:: ConjugacyClass

Cosets
++++++

.. autoclass:: Cosets

Coset
+++++

.. autoclass:: Coset

CosetTable
++++++++++

.. autoclass:: CosetTable

LeftCosets
++++++++++

.. autoclass:: LeftCosets

RightCosets
+++++++++++

.. autoclass:: RightCosets

NormalCosets
++++++++++++

.. autoclass:: NormalCosets

LeftCoset
+++++++++

.. autoclass:: LeftCoset

RightCoset
++++++++++

.. autoclass:: RightCoset

NormalCoset
+++++++++++

.. autoclass:: NormalCoset

DoubleCosets
++++++++++++

.. autoclass:: DoubleCosets

DoubleCoset
+++++++++++

.. autoclass:: DoubleCoset
