RepLAB group hierarchy
======================

.. _MonoidsGroups

Monoids and groups
++++++++++++++++++

Monoids are sets equipped with an associative binary operation and an identity; group elements have inverses.

.. autoclass:: +replab.Monoid
	       :members:
	       :show-inheritance:

.. autoclass:: +replab.Group
	       :members:
	       :show-inheritance:

.. _FinitelyGeneratedGroups

Finitely generated groups
+++++++++++++++++++++++++

All elements of finitely generated groups can be obtained by composition of a finite number of the group generators; and the set of those generators is also finite.
In particular, such groups enable the factorization of their elements into words.

.. autoclass:: +replab.FinitelyGeneratedGroup
	       :members:
	       :show-inheritance:

.. autoclass:: +replab.FreeGroup
	       :members:
	       :show-inheritance:

.. autoclass:: +replab.Word
	       :members:
	       :show-inheritance:

.. _FiniteGroups

Finite groups
+++++++++++++

Finite groups contain a finite number of elements.

.. autoclass:: +replab.FiniteGroup
	       :members:
	       :show-inheritance:

.. autoclass:: +replab.FiniteGroupDecomposition
	       :members:
	       :show-inheritance:
