Base object types and functions
===============================

The following classes are at the base of many RepLAB objects/operations.

.. module:: +replab

* `.Str`: At the root of our object hierarchy, the base class `+replab.Str` is in charge of pretty printing all the objects defined
  in RepLAB. This is important for several reasons: Matlab has some support to pretty print objects, support which is missing
  in Octave. For standard objects (matrices, vectors), there are discrepancies between the output of Matlab and Octave. Finally
  the standard display of matrices and vectors can be a bit bulky.

  We define three styles of printing for objects: a ``header`` style that prints only the size and the type of the object,
  which can be obtained using `+replab.headerStr`, a ``short`` style that fits on a single display line, which can be obtained
  using `+replab.shortStr`, and a a ``long`` style that can use multiple lines, which can be obtained using `+replab.longStr`.

* `.Domain`: describes a set of elements that can be compared for equality, and sampled from, obeying `.DomainLaws`

* `.Laws`: Most algebraic structures in RepLAB obey algebraic laws; for example the binary operation of a `.Monoid` is associative;
  and the identity element leaves other monoid elements invariant and so on. RepLAB uses randomized tests to check those axioms in
  all implementations, in the spirit of `QuickCheck <https://hackage.haskell.org/package/QuickCheck>`_
  and `ScalaCheck <https://www.scalacheck.org/>`_.
  The infrastructure is provided by the base class `.Laws`, which is returned by the `+replab.Domain.laws` method.

* `.IndexedFamily`: describes a lazily evaluated sequence of elements, with 1-based big integer indices, obeying `.IndexedFamilyLaws`

* `.dispatch`: Provides a flexible dispatch mechanism that enables the use of plugins.

Str
+++

.. autoclass:: Str

headerStr
+++++++++

.. autofunction:: headerStr

longStr
+++++++

.. autofunction:: longStr

shortStr
++++++++

.. autofunction:: shortStr

Domain
++++++

.. autoclass:: Domain

IndexedFamily
+++++++++++++

.. autoclass:: IndexedFamily

Laws
++++

.. autoclass:: Laws

dispatch
++++++++

.. autofunction:: dispatch
