Code organization in RepLAB
===========================

Directories
-----------

-  ``docs``: compiled documentation ready to be served by Jekyll,
-  ``docs_src``: source files for the executable documentation, to be
   compiled by ``replab_compiledoc.m``,
-  ``external``: external libraries, either included as a Git submodule
   (if available on a public repository such as Github), or vendored
   (like ``vpi``),
-  ``src/+replab``: the main RepLAB package, containing the main
   classes/methods,
-  ``src/+replab/+subpackages``: implementation files for different
   submodules, or experimental stuff,
-  ``tests``: tests written using MOxUnit, augmented with our laws test
   framework.

The ``Class``, ``ClassFun`` and ``ClassLaws`` triad
---------------------------------------------------

This works for ``Class`` = ``Domain``, ``Class`` = ``Semigroup``,
``Class`` = ``Action``, and so on.

There, ``Class`` is an abstract base class, while ``ClassFun`` requires
function handles that implement the abstract methods of ``Class`` as
parameters to its constructor. This enables the creation of "anonymous"
classes in the spirit of MATLAB's "anonymous" functions.

The ``ClassLaws`` class lists the algebraic laws that should be obeyed
by a particular structure.

-  ``Laws``: abstract base class for the laws checking framework.

Structured sets
---------------

Most classes in RepLAB describe sets of elements with additional
structure.

Domains
~~~~~~~

At the base of the hierarchy, the ``Domain`` describes a set of elements
that can be tested for equality (``eqv``) and from which random samples
can be taken (``sample``). Such sets are potentially infinite.

As ``Domain`` is an abstract base class, it contains abstract methods
(``eqv``, ``sample``).

To quickly create an instance of ``Domain``, the helper class
``DomainFun`` can be used, passing the method implementations as
function handles to the ``DomainFun`` constructor.

Blackbox groups
~~~~~~~~~~~~~~~

The following three classes provide basic group structures.

-  ``Monoid``: structure with a binary associatvie operation, an
   identity element (property), and test for identity
-  ``Group``: + inverse elements

Finite groups
~~~~~~~~~~~~~

Finite groups contain a finite number of elements.

-  ``FiniteGroup``: abstract base class for finite groups,
-  ``NiceFiniteGroup``: describes a finite group that has an injective
   morphism to a permutation group,

-  ``FiniteGroupDecomposition``: describes the decomposition of a group
   into a product of sets.

Abstract group constructions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``replab.directproduct.*``: represents the (outer) direct product of
   finite groups, the elements are represented by row cell arrays,
-  ``replab.semidirectproduct.*``: represents an outer semidirect
   product from a homomorphism, with both factors finite groups,
-  ``replab.wreathproduct.*``: describes the wreath product of a finite
   group by a permutation group.

Permutation groups
~~~~~~~~~~~~~~~~~~

Permutation groups are finite groups that are subgroups of the symmetric
group acting on ``n`` elements.

In particular, the decomposition of :math:`\{1,...,n\}` into orbits, the
natural action and representation of the group are all defined.

-  ``PermutationGroup``: abstract base class for permutation groups,
-  ``Permutations``: describes the symmetric group acting on ``n``
   points.

As a closely related variant, we describe signed permutations.

-  ``signed.PermutationGroup``: provides an abstract base class for all
   signed permutation groups,
-  ``signed.Permutations``: describes the signed permutations on
   :math:`\{-n,...,-1, 1,...,n\}`,

Actions
~~~~~~~

The following two classes in RepLAB do not define a domain, but rather
express a relation between domains.

-  ``Action``: abstract base class that defines the action of elements
   of a group ``G`` on elements of a domain ``P``,

Representations
~~~~~~~~~~~~~~~

-  ``Rep`` is a real or complex representation of a finitely generated
   group.

-  ``RepByImages`` is a real or complex representation described by the
   (matrix) images of the group generators; the computation of images
   rests on the factorization of group elements.

Once a ``Rep`` has been constructed, it can be decomposed into
irreducible representations using the ``realRep.decomposition`` method.

-  ``SubRep`` describes a subrepresentation of an existing
   representation. The instance can remember which representation it
   splits (as a ``parent``) and the change of basis matrix (``U``).
-  ``Irrep`` describes an irreducible subrepresentation; for
   representations over the reals, it also takes care to express the
   complex/quaternion division algebras in a regular manner.
-  ``Isotypic`` are isotypic components, which group equivalent
   irreducible representations present in a reprsentation.
-  ``Irreducible`` regroups isotypic components.

Commutant algebras
~~~~~~~~~~~~~~~~~~

The key to the simplification of invariant semidefinite programs is the
decomposition of *matrices that commute with a group representation*;
the set of all such matrices is the *commutant algebra*.

RepLAB has several classes to work with centralizer algebras.

-  ``Equivariant`` describes the vector space of equivariant linear maps
   between two representations of the same group.
-  ``Commutant`` describes a generic commutant algebra.
-  ``DivisionAlgebra`` describes a division algebra present in an
   irreducible representation.

Sequences of elements
---------------------

-  ``IndexedFamily``: describes an indexed family, whose index is a big
   integer (``vpi``).

Utility classes
---------------

-  ``Partition``: describes an unordered partition of the integers
   :math:`\{1,...,n\}`,
-  ``DivisionAlgebra``: generic representation of the real, complex and
   quaternionic algebras over the reals,
-  ``rational``: a hacky implementation of rational matrices with
   double-representable integer coefficients,
-  ``quaternion.H``: a minimal implementation of quaternion matrices,
-  ``Parameters``: various global settings such as tolerances,
-  ``GroupMorphismLaws``: laws for function handles that are group
   homomorphisms.

Methods
-------

-  ``isNonZeroMatrix``, tests whether a matrix is nonzero up to a given
   tolerance using the 2-norm (singular value), with accelerations
   provided by matrix norm inequalities and cheap norms,
-  ``headerStr``: returns a short description of an object without
   digging into its contents,
-  ``shortStr``: returns a short description of an object, with a best
   effort to include its contents,
-  ``longStr``: returns a long description of an object, including an
   enumeration of its properties.
