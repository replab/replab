Base object types in RepLAB
===========================

The following classes are at the base of many RepLAB objects.

.. module:: +replab

* `.Str`: provides infrastructure for pretty printing of objects, in a way compatible with Matlab and Octave
* `.Samplable`: describes a set of elements that can be sampled from
* `.Domain`: describes a set of elements that can be compared for equality, and sampled from

Str
+++

.. autoclass:: Str

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
