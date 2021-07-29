Abstract groups and atlas
=========================

Abstract groups are defined using a presentation, without direct reference to a concrete realization. They
are thus used when user-defined groups are recognized.

* `.AbstractGroup` is a finite group whose elements are written as words composed of products of generators.

* `.Atlas` is the infrastructure that recognizes groups. In particular, it is used to tip the user of the
  type of the group they just constructed. It is also used in exact decompositions, to retrieve the character
  table of a particular finite group.

* `.AtlasEntry` contains a sketch of the structure of a finite group known by RepLAB. It is used to speed up
  the recognition process.

.. module:: +replab

AbstractGroup
+++++++++++++

.. autoclass:: AbstractGroup

Atlas
+++++

.. autoclass:: Atlas

AtlasEntry
++++++++++

.. autoclass:: AtlasEntry
