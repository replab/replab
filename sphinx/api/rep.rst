Representations
===============

.. module:: +replab

-  `.Rep` is a real or complex representation of a finitely generated
   group with associated `.RepLaws`.

-  `.RepByImages` is a real or complex representation described by the
   (matrix) images of the group generators; the computation of images
   rests on the factorization of group elements.

-  `.SimilarRep` describes a representation with an associated change of basis matrix.

Once a `.Rep` has been constructed, it can be decomposed into
irreducible representations using the ``rep.decomposition`` method.

-  `.SubRep` describes a subrepresentation of an existing
   representation. The instance can remember which representation it
   splits (as a ``parent``) and the change of basis matrix (``basis``).

-  `.Isotypic` are isotypic components, which group equivalent
   irreducible representations present in a representation. If the irreducible
   representations are all in the same basis, this will be an instance of
   `.HarmonizedIsotypic`.

-  `.Irreducible` regroups isotypic components, with associated
   `.IrreducibleLaws`.

Representations induce vector spaces.

-  `.Equivariant` describes the vector space of equivariant linear maps
   between two representations of the same group, with specialization
   for commutant maps of isotypic and irreducible decompositions:
   `.IsotypicCommutant` and `.IrreducibleCommutant`.

Rep
+++

.. autoclass:: Rep

RepByImages
+++++++++++

.. autoclass:: RepByImages

SimilarRep
++++++++++

.. autoclass:: SimilarRep

SubRep
++++++

.. autoclass:: SubRep

Isotypic
++++++++

.. autoclass:: Isotypic

HarmonizedIsotypic
++++++++++++++++++

.. autoclass:: HarmonizedIsotypic

Irreducible
+++++++++++

.. autoclass:: Irreducible

Equivariant
+++++++++++

.. autoclass:: Equivariant

IsotypicCommutant
+++++++++++++++++

.. autoclass:: IsotypicCommutant

IrreducibleCommutant
++++++++++++++++++++

.. autoclass:: IrreducibleCommutant

CharacterTable
++++++++++++++

.. autoclass:: CharacterTable
