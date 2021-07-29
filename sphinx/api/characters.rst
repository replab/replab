Characters
==========

.. module:: +replab

Character table based methods are used in RepLAB to provide exact results when requested, at the price of a major slowdown.

- `.Character` represents a character of a finite group (not necessarily irreducible).

- `.CharacterTable` is the base class

- `.ComplexCharacterTable` is the usual textbook character table of a finite group.

- `.RealCharacterTable` is the character table of the real irreducible representations of a finite group.

See also `.ConjugacyClasses` and `.ConjugacyClass`.

Character
+++++++++

.. autoclass:: Character

CharacterTable
++++++++++++++

.. autoclass:: CharacterTable

ComplexCharacterTable
+++++++++++++++++++++

.. autoclass:: ComplexCharacterTable

RealCharacterTable
++++++++++++++++++

.. autoclass:: RealCharacterTable
