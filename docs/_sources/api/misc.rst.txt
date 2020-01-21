Miscellaneous functions/classes
===============================

.. module:: +replab

-  `.Partition`: describes an unordered partition of the integers :math:`\{1,...,n\}`,

-  ``settings``: various global settings such as tolerances,

-  `.isNonZeroMatrix`, tests whether a matrix is nonzero up to a given
   tolerance using the 2-norm (singular value), with accelerations
   provided by matrix norm inequalities and cheap norms.

Partition
+++++++++

.. autoclass:: Partition

isNonZeroMatrix
+++++++++++++++

.. autofunction:: isNonZeroMatrix
