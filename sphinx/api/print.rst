Pretty printing
===============

At the root of our object hierarchy, the base class `+replab.Str` is in charge of pretty printing all the objects defined
in RepLAB. This is important for several reasons: Matlab has some support to pretty print objects, support which is missing
in Octave. For standard objects (matrices, vectors), there are discrepancies between the output of Matlab and Octave. Finally
the standard display of matrices and vectors can be a bit bulky.

We define three styles of printing for objects:

* a ``header`` style that prints only the size and the type of the object,
* a ``short`` style that fits on a single display line,
* a ``long`` style that can use multiple lines.

The style used when displaying objects in the REPL/command line is the ``long`` style.

The following classes and functions are at the base of the RepLAB pretty printing framework.

.. module:: +replab

* `.Str`: provides infrastructure for pretty printing of objects, in a way compatible with Matlab and Octave
* `+replab.headerStr`, `+replab.longStr`, `+replab.shortStr`: Prints arbitrary objects
* ... the rest of the elements below are helpers in the ``replab.str`` package


Str
+++

.. autoclass:: Str

headerStr
+++++++++

.. autofunction:: headerStr

longStr
+++++++

.. autofunction:: longStr

shortStr
++++++++

.. autofunction:: shortStr

.. module:: +replab.+str

str.align
+++++++++

.. autofunction:: align

str.brackets
++++++++++++

.. autofunction:: brackets

str.cellStr
+++++++++++

.. autofunction:: cellStr

str.escape
++++++++++

.. autofunction:: escape

str.fieldsList
++++++++++++++

.. autofunction:: fieldsList

str.headerStr
+++++++++++++

.. autofunction:: headerStr

str.horzcatForce
++++++++++++++++

.. autofunction:: horzcatForce

str.longFit
+++++++++++

.. autofunction:: longFit

str.longStr
+++++++++++

.. autofunction:: longStr

str.Normalized
++++++++++++++

.. autoclass:: Normalized

str.pluralize
+++++++++++++

.. autofunction:: pluralize

str.shortStr
++++++++++++

.. autofunction:: shortStr

str.sizeStr
+++++++++++

.. autofunction:: sizeStr

str.uniqueNames
+++++++++++++++

.. autofunction:: uniqueNames
