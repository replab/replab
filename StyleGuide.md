# Reuse of primitive types

- Permutations are 1 x n vectors containing a permutation of the integers {1..n}, stored as double coefficients.
- Signed permutations are 1 x n vectors containing the integers {+/-1..+/-n}, such that `abs(signedPerm)` is a permutation.

# Class structure

For each algebraic structure (semigroup, monoid, group, ...), we define:

- An abstract base class (Semigroup) with no constructor, but whose
  abstract methods try to call handle function present as properties;
  this is used in conjunction with the XXXFun implementation
- A law checking class (SemigroupLaws)
- A function handle quick implementation class SemigroupFun
  that has a single constructor, calling all XXXFun superclasses

# Abstract classes/methods

As Octave does not support abstract methods, we use the following workaround.

In the base abstract class, an abstract method `operation(args)` looks for a property called `operationFun`, stores it into a variable `f` then calls `f(args)`.

This enables subclasses to implement the abstract class in a generic way (see `replab.GroupFun` for an example) without duplicating too much code.

If an error appears at runtime in one of those pseudo-abstract methods, it was supposed to be overriden in a subclass but was not.

# Typeclasses

Whenever possible, we reuse primitive types and do not wrap them in classes.

To perform generic operations upon them, we pass around a "typeclass", which is a collection of methods/default values for the type.

All typeclasses derive from `replab.Laws` which provides the test harness for [property-based checks](https://en.wikipedia.org/wiki/QuickCheck).

## Domains

Some typeclasses define *domains*, which are mathematical sets.
Domains define the following methods:

- `eqv` tests for equality,
- `sample` provides a random element of the domain.

Example: `replab.Permutations(10)` is the domain of permutations acting on 10 elements.

## Non-domain typeclasses

Other typeclasses define relations between domains.

## Type properties

Typeclasses that depend on other types contain the typeclasses of those types as properties.
For example, an `Action` typeclass contain a property `G` describing the group structure, and a property `P` describing the domain of elements being acted upon.
There, `G` and `P` are domains, while `Action` itself is not.

## Laws

Typeclasses can define laws, which are defined as methods in a companion class deriving from `replab.Laws`. These methods are of the form `law_name_of_the_law_TYPES`, where `TYPES` describes the parameters of the law (see `replab.Laws`).
