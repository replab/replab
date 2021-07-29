Equivariant spaces
==================

.. module:: +replab

Equivariant spaces describe linear map between representations that preserve the group structure.

The hierarchy of equivariant spaces closely mirrors the hierarchy of representations `.Rep`.

- `.Equivariant` describes the vector space of equivariant linear maps between two representations of the same group.

- `.SubEquivariant` describes a subspace of an equivariant space induced by subspaces of source/target representations.

- `.IsotypicEquivariant` describes a subspace of an equivariant space corresponding to a pair of isotypic subspaces of the source/target representations

- `.IrreducibleEquivariant` is the most interesting class, as it describes the equivariant space in its block form, taking in account the decomposition of the source and target representation into irreps.

Equivariant
+++++++++++

.. autoclass:: Equivariant

SubEquivariant
++++++++++++++

.. autoclass:: SubEquivariant

IsotypicEquivariant
+++++++++++++++++++

.. autoclass:: IsotypicEquivariant

IrreducibleEquivariant
++++++++++++++++++++++

.. autoclass:: IrreducibleEquivariant
