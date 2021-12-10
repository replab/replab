Group morphisms
===============

.. module:: +replab

Group morphisms are used in **RepLAB** to relate groups to one another. For example, there is always a (trivial) injection map that embeds a subgroup into the group that contains it.

A morphism is simply a function between groups that preserve the multiplicative structure ($f(x y) = f(x) f(y)$). An isomorphism is an invertible (or bijective) morphism.

Isomorphisms are used in **RepLAB** to compute with finite groups, by translating the operations into an isomorphic permutation group.

Isomorphisms are also used when the user asks **RepLAB** to recognize a finite group (`.FiniteGroup.recognize`). **RepLAB** then returns the isomorphism from the abstract group to the user-defined group. Thus, the generators of the "platonic" abstract group can be identified as elements of the user-defined group.

Morphisms are used to describe how symmetry groups relate to one another in the convex optimization framework (see `.equiop`).

* `.Morphism` is the most general class and a base for all the classes below.

* `.FiniteMorphism` is a morphism between finite groups; in particular, it allows the computation of preimages.

* `.Isomorphism` is a bijective morphism.

* `.FiniteIsomorphism` is a bijective finite morphism.

Morphism
++++++++

.. autoclass:: Morphism

Isomorphism
+++++++++++

.. autoclass:: Isomorphism

FiniteMorphism
++++++++++++++

.. autoclass:: FiniteMorphism

FiniteIsomorphism
+++++++++++++++++

.. autoclass:: FiniteIsomorphism
