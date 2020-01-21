Style guide
===========

We describe guiding principles used in the code.

Capitalization
--------------

CamelCase/camelCase follows the separations of US English; words
separated either by a hyphen or a space are separated in the camel case.
Class names are ``UpperCamelCase``. Method and function names are
``lowerCamelCase``.

Object creation
---------------

When there is only one implementation of an interface, the constructor
is public and can be called directly.

Example: ``S4 = replab.Permutations(4)``

When a base interface is available and several optimized implementations
are available, the object should be constructed by calling a function
with the same name as the base interface, in ``lowerCamelCase``.

Example: ``c = replab.equivariant(rep1, rep2)``

That function will then construct the object with the best performance
characteristics.

No abstract methods in abstract base classes
--------------------------------------------

As Octave does not support the ``methods (Abstract)`` syntax, we provide
generic implementations for abstract methods by throwing
``error('Abstract')``.
