Base group structures
=====================

.. module:: +replab

The following classes provide basic structures for monoids/groups. All objects below inherit from `.Str` / `.Obj` / `.Domain`, and
thus provide pretty printing infrastructure, a `.Domain.eqv` equality test method, and a `.Domain.sample` method to generate a random element.

* `.Monoid`: Structure with a binary associative operation, an identity element (property), and test for identity.
* `.Group`: A monoid with inverse elements.
* `.CompactGroup`: A `.Group` with an associated Haar measure, which enables the construction of representations that can be handled by RepLAB.
* `.AutomorphismGroup`: Describes the automorphism group of a group
* `.Action`: Group can act on objects, and this describes the action of elements of a group on elements of a domain.

Monoid
++++++

.. autoclass:: Monoid

Group
+++++

.. autoclass:: Group

CompactGroup
++++++++++++

.. autoclass:: CompactGroup

AutomorphismGroup
+++++++++++++++++

.. autoclass:: AutomorphismGroup

Action
++++++

.. autoclass:: Action
