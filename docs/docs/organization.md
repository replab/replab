# Code organization in RepLAB

## Directories

- `docs`: compiled documentation ready to be served by Jekyll,
- `docs_src`: source files for the executable documentation, to be compiled by `replab_compiledoc.m`,
- `external`: external libraries, either included as a Git submodule (if available on a public repository such as Github), or vendored (like `vpi`),
- `+replab`: the main RepLAB package, containing the main classes/methods,
- `+replab/+subpackages`: implementation files for different submodules,
- `tests`: tests written using MOxUnit, augmented with our laws test framework.

## No abstract methods in abstract base classes

As Octave does not support the `methods (Abstract)` syntax, we provide generic implementations for abstract methods in abstract base classes by either:

- trying to retrieve and calling a function handle property (the `DomainFun`, `SemigroupFun`, ... approach), or
- throwing `error('Not implemented')`.

## The `Class`, `ClassFun` and `ClassLaws` triad

This works for `Class` = `Domain`, `Class` = `Semigroup`, `Class` = `Action`, and so on.

There, `Class` is an abstract base class, while `ClassFun` requires function handles that implement the abstract methods of `Class` as parameters to its constructor. This enables the creation of "anonymous" classes in the spirit of MATLAB's "anonymous" functions.

The `ClassLaws` class lists the algebraic laws that should be obeyed by a particular structure.

- `Laws`: abstract base class for the laws checking framework.

## Pretty printing

Most all classes in RepLAB derive from the `Str` base class, which implements a `str` method that provides a text representation of a RepLAB object.

RepLAB provides also a `strOf` function that calls the `str` method for instances of `Str`, and implements reasonable defaults for objects outside RepLAB (such as matrices / vectors).

The `Str` class provides a default `disp` method.

The default `str` method returns the content of the `description` property (if nonempty), or constructs a text representation.

## Structured sets

Most classes in RepLAB describe sets of elements with additional structure.

### Domains

At the base of the hierarchy, the `Domain` describes a set of elements that can be tested for equality (`eqv`) and from which random samples can be taken (`sample`). Such sets are potentially infinite.

As `Domain` is an abstract base class, it contains abstract methods (`eqv`, `sample`).

To quickly create an instance of `Domain`, the helper class `DomainFun` can be used, passing the method implementations as function handles to the `DomainFun` constructor.

### Blackbox groups

The following three classes provide basic group structures.

- `Semigroup`: structure with associative binary operation
- `Monoid`: + an identity element (property), and test for identity
- `Group`: + inverse elements

### Finitely generated groups

All elements of finitely generated groups can be obtained by composition of a finite number of the group generators; and the set of those generators is also finite.

In particular, such groups enable the factorization of their elements into words.

- `FinitelyGeneratedGroup`: abstract base class for finitely generated groups,
- `FreeGroup`: the free group on `n` generators,
- `Word`: class describing elements of `FreeGroup`.

### Finite groups

Finite groups contain a finite number of elements.

- `FiniteGroup`: abstract base class for finite groups,
- `SignedPermutations`: describes the signed permutations on $\{-n,...,-1, 1,...,n\}$.

- `FiniteGroupDecomposition`: describes the decomposition of a group into a product of sets.

### Permutation groups

Permutation groups are finite groups that are subgroups of the symmetric group acting on `n` elements.

In particular, the decomposition into $\{1,...,n\}$ into orbits, the natural action and representation of the group are all defined.

- `PermutationGroup`: abstract base class for permutation groups,
- `Permutations`: describes the symmetric group acting on `n` points.

### Actions

Two classes in RepLAB do not define a domain, but rather express a relation between two domains.

- `Action`: abstract base class that defines the action of elements of a group `G` on elements of a domain `P`,
- `FaithfulAction`: in addition, for any `g` in `G` which is not the identity, the class provides a `p` in `P` that is moved by `g`.

### BSGS groups

Base and strong generating set structures enable powerful algorithms to study finite groups.

Those algorithms are generic and rest on the availability of a `FaithfulAction` for the group.

- `BSGSGroup`: generic implementation,
- `perm.PermutationBSGSGroup`: implementation for permutation groups.

## Sequences of elements

- `Enumerator`: describes an indexed sequence, whose index is a big integer (`vpi`).

## Utility classes

- `DivisionAlgebra`: generic representation of the real, complex and quaternionic algebras over the reals,
- `rational`: a hacky implementation of rational matrices with double-representable integer coefficients,
- `Settings`: various global settings such as tolerances,
- `GroupMorphismLaws`: laws for function handles that are group homomorphisms.

## Methods

- `isNonZeroMatrix`, tests whether a matrix is nonzero up to a given tolerance using the 2-norm (singular value), with accelerations provided by matrix norm inequalities and cheap norms,
- `strOf`: prettyprints a generic object.
