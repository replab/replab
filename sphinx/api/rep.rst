Representations
===============

.. module:: +replab

Representation in RepLAB, whether they are exact or inexact, reducible or irreducible, are described by the `.Rep` base class.

* `.Rep` is a real or complex representation of a compact group.

Some representations encode a division algebra using real or complex coefficients.

* `.DivisionAlgebra` provides information about the division algebra encoding handled by RepLAB.

A particular type of representation is a representation of a finite group given by generator images.

*  `.RepByImages` is a real or complex representation described by the (matrix) images of the group generators.

Once a `.Rep` has been constructed, it can be decomposed into irreducible representations using the ``rep.decomposition`` method.

*  `.SubRep` describes a subrepresentation of an existing representation. The instance can remember which
   representation it splits (`.SubRep.parent`) and the change of basis maps (`.SubRep.injection`, `.SubRep.projection`).

*  `.Isotypic` are isotypic components, which group equivalent irreps present in a representation.

*  `.Irreducible` represents the decomposition of a representation into isotypic components.

Rep
+++

.. autoclass:: Rep

DivisionAlgebra
+++++++++++++++

.. autoclass:: DivisionAlgebra

RepByImages
+++++++++++

.. autoclass:: RepByImages

SubRep
++++++

.. autoclass:: SubRep

Isotypic
++++++++

.. autoclass:: Isotypic

Irreducible
+++++++++++

.. autoclass:: Irreducible
