RepLAB pretty printing infrastructure
=====================================

Matlab has little support for consistent pretty printing of objects.
Thus, RepLAB contains helper classes and functions to print arbitrary objects in
a clear style. This :class:`+replab.Str` base class provides sane defaults
and is used by all RepLAB classes. 

Three styles of printing are used in RepLAB.

* a *header* style that prints only the size and type of the object,

* a *short* style that fits on a single display line,

* a *long* style that can use multiple lines.

The style used when displaying objects in the REPL/command line is the
*long* style. When displaying the properties embedded of an object,
the *short* style is used. The *header* style is used when the short style
does not fit in the maximal number of column characters.
  
When used on a :class:`+replab.Str` instance, the RepLAB functions
:func:`+replab.headerStr`, :func:`+replab.shortStr`, :func:`+replab.longStr`
use the methods of this base class. Additionally, they provide sane defaults when
called on other Matlab types, and have default values for the maximal number of
displayed rows/columns.

.. _StrMethods:

Pretty printing methods
+++++++++++++++++++++++

Depending on the desired verbosity, you can call either of :func:`+replab.headerStr`,
:func:`+replab.shortStr` or :func:`+replab.longStr` to obtain the string
representation of an arbitrary Matlab object.

.. autofunction:: +replab.headerStr

.. autofunction:: +replab.shortStr

.. autofunction:: +replab.longStr

.. _Str:

Str base class
++++++++++++++

.. autoclass:: +replab.Str
     :members:
     :show-inheritance:

.. _StrHelpers:

Pretty printing helper methods
++++++++++++++++++++++++++++++

.. autofunction:: +replab.+str.align
