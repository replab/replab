Base group structures
=====================

.. module:: +replab

The following classes provide basic structures for monoids/groups. All objects below inherit from `.Str` and `.Domain`, and
thus provide pretty printing infrastructure, a `.Domain.eqv` equality test method, and a `.Domain.sample` method to generate a random element.

* `.Monoid`: Structure with a binary associative operation, an identity element (property), and test for identity.
* `.Group`: Adds inverse elements.
* `.CompactGroup`: A `.Group` with an associated Haar measure, which enables the construction of representations that can be handled by RepLAB.
* `.FiniteGroup`: A group with a finite number of elements, is generated by a finite number of generators. Moreover, it can be decomposed in a cartesian product of sets to speed up averaging (`.FiniteGroupDecomposition`).
* `.NiceFiniteGroup`: Describes a group with an injective homomorphism onto a permutation group; it is a base class of other groups.

Group can have group actions (`.Action`), which is the action of elements of a group ``G`` on elements of a domain ``P``,

Monoid
++++++

.. autoclass:: Monoid

Group
+++++

.. autoclass:: Group

CompactGroup
++++++++++++

.. autoclass:: CompactGroup

FiniteGroup
+++++++++++

.. autoclass:: FiniteGroup

FiniteGroupDecomposition
++++++++++++++++++++++++

.. autoclass:: FiniteGroupDecomposition

NiceFiniteGroup
+++++++++++++++

.. autoclass:: NiceFiniteGroup

Action
++++++

.. autoclass:: Action
