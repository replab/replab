Permutation groups
==================

.. module:: +replab

RepLAB provides ways to work with permutation groups.

Permutation groups are finite groups that are subgroups of the symmetric group acting on ``n`` elements.

In particular, the decomposition of :math:`\{1,...,n\}` into orbits, the natural action and representation of the group are all defined.

* `.Permutation`: Class with static methods to work with permutations.
* `.PermutationGroup`: Describes a generic permutation group.

As a closely related variant, we describe signed permutations. Signed permutations act on the domain
:math:`\{-n,...,-1, 1,...,n\}`.

* `.SignedPermutation`: Class with static methods to work with signed permutations.
* `.SignedPermutationGroup`: Describes a group of signed permutations.


Permutation
+++++++++++

.. autoclass:: Permutation

PermutationGroup
++++++++++++++++

.. autoclass:: PermutationGroup

SignedPermutation
+++++++++++++++++

.. autoclass:: SignedPermutation

SignedPermutationGroup
++++++++++++++++++++++

.. autoclass:: SignedPermutationGroup
