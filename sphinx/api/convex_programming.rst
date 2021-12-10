Convex programming
==================

.. module:: +replab

A major application of **RepLAB** is symmetrizing convex problems to compute their optimum faster.

In particular, it supports:

* `.CommutantVar` which provides a substitute for YALMIP's ``sdpvar``, with a lot of helper methods. `.CommutantVar` will be replaced in the future by a more comprehensive framework.

* `.SedumiData` which preprocesses SDP problems in an extension of the format used by SeDuMi

* `.equivar` and `.equiop` are the next version of our symmetric convex optimization framework. They are still experimental.

CommutantVar
++++++++++++

.. autoclass:: CommutantVar

SedumiData
++++++++++

.. autoclass:: SedumiData

equivar
+++++++

.. autoclass:: equivar

equiop
++++++

.. autoclass:: equiop
