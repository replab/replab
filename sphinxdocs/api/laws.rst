Laws
====

.. module:: +replab

Most algebraic structures in RepLAB obey algebraic laws; for example the binary operation of a `.Monoid` is associative;
and the identity element leaves other monoid elements invariant and so on.

RepLAB uses randomized tests to check those axioms in all implementations, in the spirit of `QuickCheck <https://hackage.haskell.org/package/QuickCheck>`_
and `ScalaCheck <https://www.scalacheck.org/>`_.


The infrastructure is provided by the base class `.Laws`.

.. _Laws:

Laws
++++

.. autoclass:: Laws
