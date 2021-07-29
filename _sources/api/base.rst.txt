Base object types and functions
===============================

The following classes are at the base of many RepLAB objects/operations.

.. module:: +replab

* `.Str`: Inherited by all RepLAB objects, the base class `+replab.Str` is in charge of pretty printing.
  This is important for several reasons: Matlab has some support to pretty print objects, support which is missing in Octave.
  For standard objects (matrices, vectors), there are discrepancies between the output of Matlab and Octave.
  Finally the standard display of matrices and vectors can be a bit bulky.

  We define three styles of printing for objects, which are covered by the three functions below.

* `.headerStr`: a ``header`` style that prints only the size and the type of the object,

* `.shortStr`: a ``short`` style that fits on a single display line

* `.longStr`: a ``long`` style that uses multiple lines, which can be obtained using `+replab.longStr`.

* `.Obj`: This class is at the base of most RepLAB objects with notable exceptions (`.cyclotomic` and `.H`, which are value classes).
  The `.Obj` base class provides helper methods to check the laws obeyed by the object in question (`.Obj.laws` and `.Obj.check`),
  methods to cache computed properties, and finally a default implementation of object comparison using ``==``.

* `.Laws`: Most algebraic structures in RepLAB obey algebraic laws; for example the binary operation of a `.Monoid` is associative;
  and the identity element leaves other monoid elements invariant and so on. RepLAB uses randomized tests to check those axioms in
  all implementations, in the spirit of `QuickCheck <https://hackage.haskell.org/package/QuickCheck>`_
  and `ScalaCheck <https://www.scalacheck.org/>`_.
  The infrastructure is provided by the base class `.Laws`, which is returned by the `+replab.Obj.laws` method.

* `.Domain`: describes a set of elements that can be compared for equality, and randomly sampled from. The capability of sampling random
  elements from a domain is heavily used in law checks.

Note that some objects implement standard Matlab equality semantics (a.k.a. ``x == y``); for example cyclotomic matrices can be tested
for equality, and finite groups of the same type as well. Other objects such as generic compact groups, representations, will be compared
using "is this the same object in memory".

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

Obj
+++

.. autoclass:: Obj

Domain
++++++

.. autoclass:: Domain

Laws
++++

.. autoclass:: Laws
