Style guide
===========

.. module:: +replab

We describe guiding principles used in the code.

No continuations in method/function/class declarations
------------------------------------------------------

Lines with a ``function`` or ``classdef`` declaration must be standalone, i.e. must not end with a continuation ``...``.

This mean that those lines can go over the column limit (TODO: add link to this limit when we decide on one.)

Reuse of primitive types
------------------------

Whenever possible, we reuse primitive types and do not wrap them in classes.

To perform generic operations upon them, we pass around a "typeclass", which is a collection of methods/default values for the type.

All typeclasses derive from `.Laws` which provides the test harness for `property-based checks <https://en.wikipedia.org/wiki/QuickCheck>`_.

- Permutations are ``1 x n`` vectors containing a permutation of the integers ``{1..n}``, stored as double coefficients.

- Signed permutations are ``1 x n`` vectors containing the integers ``{+/-1..+/-n}``, such that ``abs(signedPerm)`` is a permutation.

Typeclasses
-----------

Some typeclasses define *domains*, which are mathematical sets.

Domain defining typeclasses
...........................

Domains define the following methods:

- `+replab.+Domain.eqv` tests for equality,
- `+replab.+Domain.sample` provides a random element of the domain.

Example: ``replab.Permutations(10)`` is the domain of permutations acting on 10 elements.

Other typeclasses
.................

Other typeclasses define relations between domains, such as `+replab.Action`.

Typeclasses that depend on other types contain the typeclasses of those types as properties.
For example, a `+replab.Action` typeclass contain a property ``G`` describing the group structure, and a property ``P`` describing the domain of elements being acted upon.
There, ``G`` and ``P`` are domains, while `~replab.Action` itself is not.

Laws
....

Typeclasses can define laws, which are defined as methods in a companion class deriving from `+replab.Laws`. These methods are of the form ``law_name_of_the_law_TYPES``, where ``TYPES`` describes the parameters of the law (see `+replab.Laws`).

Class structure
---------------

For each algebraic structure (monoid, group, ...), we define:

- An abstract base class (ex: `.Monoid`) with no constructor, where the abstract methods have
  a body containing an error (see below).

- A law checking class (ex: `.MonoidLaws`)

- (optional) A generic implementation using function handles (ex: `+replab.+lambda.Monoid`)

Capitalization
--------------

CamelCase/camelCase follows the separations of US English; words
separated either by a hyphen or a space are separated in the camel case.
Class names are ``UpperCamelCase``. Method and function names are
``lowerCamelCase``; exception: functions that create objects such as `DihedralGroup`.

Lazy properties
---------------

If a computed property is cached, it is implemented using a protected property `expensiveProp_`, and accessed using a method
`expensiveProp` that then checks if `expensiveProp_` has already been computed.

Object creation
---------------

When there is only one implementation of an interface, the constructor
is public and can be called directly.

Example: ``S4 = replab.Permutations(4)``

When a base interface is available and several optimized implementations
are available, the object should be constructed by calling a function
with the same name as the base interface with a ``make`` prefix, in ``lowerCamelCase``.

Example: ``c = replab.makeEquivariant(rep1, rep2)``

That function will then construct the object with the best performance
characteristics.

No abstract methods in abstract base classes
--------------------------------------------

As Octave does not support the ``methods (Abstract)`` syntax, we provide
generic implementations for abstract methods by having the method body be
a single line with ``error('Abstract');``.

Use of cell arrays vs. varargin
-------------------------------

In general, we avoid using methods/functions with a variable number of arguments, as it makes future extensions with optional parameters difficult. Exception: helper methods/function names with ``of`` in their names.
