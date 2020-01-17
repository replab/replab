Base object types in RepLAB
===========================

The following classes are at the base of many RepLAB objects.

.. module:: +replab

* `.Samplable`: describes a set of elements that can be sampled from
* `.Domain`: describes a set of elements that can be compared for equality, and sampled from, obeying `.DomainLaws`
* `.IndexedFamily`: describes a lazily evaluated sequence of elements, with 1-based big integer indices, obeying `.IndexedFamilyLaws`

Samplable
+++++++++

.. autoclass:: Samplable

Domain
++++++

.. autoclass:: Domain

DomainLaws
++++++++++

.. autoclass:: DomainLaws

IndexedFamily
+++++++++++++

.. autoclass:: IndexedFamily

IndexedFamilyLaws
+++++++++++++++++

.. autoclass:: IndexedFamilyLaws
