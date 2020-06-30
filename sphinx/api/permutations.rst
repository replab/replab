Permutation groups
==================

.. module:: +replab

RepLAB provides ways to work with permutation groups.

Permutation groups are finite groups that are subgroups of the symmetric
group acting on ``n`` elements.

In particular, the decomposition of :math:`\{1,...,n\}` into orbits, the
natural action and representation of the group are all defined.

* `.PermutationGroup`: provides the infrastructure for all permutation groups, with associated `.PermutationGroupLaws`,
* `.Permutations`: describes the symmetric group on a given domain, with associated `.PermutationsLaws`,
* `.S`: Shortcut for `.Permutations`.

As a closely related variant, we describe signed permutations. Signed permutations act on the domain
:math:`\{-n,...,-1, 1,...,n\}`.

*  `.SignedPermutationGroup`: provides a base class for all signed permutation groups,
*  `.SignedPermutations`: describes the signed permutation group, with associated `signed.PermutationsLaws`.


PermutationGroup
++++++++++++++++

.. autoclass:: PermutationGroup

PermutationGroupLaws
++++++++++++++++++++

.. autoclass:: PermutationGroupLaws


Permutations
++++++++++++

.. autoclass:: Permutations

S
+

.. autofunction:: S

PermutationsLaws
++++++++++++++++

.. autoclass:: PermutationsLaws

SignedPermutationGroup
++++++++++++++++++++++

.. autoclass:: SignedPermutationGroup

SignedPermutations
+++++++++++++++++++

.. autoclass:: Permutations

SignedPermutationsLaws
+++++++++++++++++++++++

.. autoclass:: PermutationsLaws
