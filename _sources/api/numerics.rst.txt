Numerics
========

.. module:: +replab

RepLAB provides the following multi-dimensional array types as an extension to the types provided by MATLAB/Octave.

* `.cyclotomic` implements exact computation on the cyclotomic field. This enables RepLAB to handle any irreducible representation of a finite group, if care is taken to express it over the cyclotomics (which is always possible).

* `.H` implements matrices of quaternions using a pair of complex matrices. Quaternion matrices are used to decompose real representations of quaternion-type.

cyclotomic
++++++++++

.. autoclass:: cyclotomic

H
+

.. autoclass:: H
