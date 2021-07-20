Representations
===============

.. module:: +replab

-  `.Rep` is a real or complex representation of a finitely generated
   group with associated `.RepLaws`.

-  `.RepByImages` is a real or complex representation described by the
   (matrix) images of the group generators; the computation of images
   rests on the factorization of group elements.

Once a `.Rep` has been constructed, it can be decomposed into
irreducible representations using the ``rep.decomposition`` method.

-  `.SubRep` describes a subrepresentation of an existing
   representation. The instance can remember which representation it
   splits (as a ``parent``) and the change of basis matrix (``basis``).

-  `.Isotypic` are isotypic components, which group equivalent
   irreducible representations present in a representation.

-  `.Irreducible` regroups isotypic components, with associated
   `.IrreducibleLaws`.

Representations induce vector spaces.

-  `.Equivariant` describes the vector space of equivariant linear maps
   between two representations of the same group.

Rep
+++

.. autoclass:: Rep

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

Equivariant
+++++++++++

.. autoclass:: Equivariant

CharacterTable
++++++++++++++

.. autoclass:: CharacterTable
