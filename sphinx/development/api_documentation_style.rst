API Documentation style guide
=============================

This document describes the syntax and best practices for the
documentation comments embedded in the source code.

Overview
--------

This document is heavily inspired by the `numpydoc docstring
guide <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`__
and follows its structure.

We employ the variants proposed by the `Google
style <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`__
for compactness.

We aim to achieve the following middle ground:

-  The documentation comments can be processed easily by the embedded
   ``help`` function to present very readable and complete documentation
   to the matlab prompt.

-  The documentation comments can be parsed with Sphinx to produce
   stylish and naviguable documentation.

Documentation comments
----------------------

MATLAB documentation calls the documentation strings that appears at the
beginning of a class/function/method/property definition ``comments``.
To distinguish them from code comments, we call those special comments
``documentation comments``, noting that the ``docstring`` terminology is
specific to the Python language (standalone strings in the source code).

The documentation comments are written using a subset of the
`reStructuredText <http://docutils.sourceforge.net/rst.html>`__ syntax.
This style guide does not list which syntax elements are permitted: use
common sense to employ the conventions that are compatible with
reStructuredText, and display well standalone on the terminal (as is
advocated by the numpydoc guide linked above).

In particular, we use single backticks for class, function, method,
property names, and double backticks for verbatim text. Moreover, we
include the ``+`` symbol in from of packages.

    A guiding principle is that human readers of the text are given
    precedence over contorting docstrings so our tools produce nice
    output.

The length of comment lines should be kept to 75 characters; note that
this applies to the comment string itself, not including the ``%``
character, whitespace, and source code tokens. RepLAB uses longer source
code lines.

First part of a documentation comment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a class, function or method (uniquely referred to as function from
now on), a documentation comment starts with the following introductory
content:

1. First line with a concise description of the function. Names are not
   capitalized.

2. Second line empty

3. Third line starts an extended description (and/or discussion) of the
   function follows the third line. This should clarify the role played
   by the function, i.e. its *functionality*, not to discuss
   implementation details or background theory, which should rather be
   explored in the **Notes** section below. You may refer to the
   parameters and the function name (using single backticks), but the
   parameter descriptions belong in the **Parameters** section.

Here is an example of these first few lines of documentation comment:

.. code:: matlab

    function suite = fromMethod(testClass, testMethod, varargin)
    % Create a suite from a single test method
    %
    % Creates a TestSuite from the test class described by testclass and the
    % test method described by testmethod and returns it in suite. testclass
    % is a meta.class instance which describes the desired test class. The
    % ...

When the concise description is sufficiently clear to define the
function, the complete description can be omitted. Here is an example:

.. code:: matlab

    function z = plus(x, y)
    % Addition operation

When no one-line description of the function is present, the first line
of the extended description will be used as a summary description
(possibly incomprehensible).

After these initial informations come the sections described below,
providing further systematic details.

Sections
~~~~~~~~

The documentation comment can contain a number of sections separated by
headings. Each heading is given by a keyword from the following list,
with the section content indented.

We use the following section headings:

-  ``Args`` as an alias to Sphinx ``Parameters``

-  ``Returns`` (for multiple return arguments, see the use of underlined
   section headers below).

-  ``Raises``

-  ``Example``

-  ``See Also``

-  ``Notes``

See also the
`Napoleon <https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html>`__
documentation for other possible section headers.

In order, the sections of a function comment are:

1. Arguments/Parameters

   Description of the function arguments, keywords and their respective
   types, using Google style. Argument names are not capitalized.

   .. code:: matlab

      % Args:
      %   x (type): Description of parameter ``x``.
      %   y: Description of parameter ``y`` (with type not specified)

   If it is not necessary to specify an argument, use ``optional``:

   .. code:: matlab

      %   x (type, optional): optional parameter of type ``type``

   When a parameter can only assume one of a fixed set of values, those
   values can be listed in braces. If the argument is optional, the default
   value appears first:

   .. code:: matlab

      %   field ({'R', 'C'}): Field name

   When two or more input parameters have exactly the same type, shape and
   description, they can be combined:

   .. code:: matlab

      %   a, b (double): Elements to sum

   Matlab/Octave encodes integer values using double floating point
   numbers. Our convention is to name those values ``integer``. In the rare
   case a primitive integer type is needed, we write that type precisely
   (as in ``int32`` or ``uint32``). Big integers have type ``vpi``, which
   is the name of the external library that supports them.

   Permutations are stored using row double vectors containing integers,
   and are documented as ``permutation``. Same convention for signed
   permutations, that are documented as ``signed permutation``.

   Strings represented as char arrays have type ``charstring`` (recent
   Matlab versions have a new ``string`` type which should not be confused
   with).

   Function handles have type ``function_handle``.

2. Returns

   Explanation of the returned values and their types. Return values are
   not capitalized. We distinguish two cases.

   First, when a single value is returned, we use the Google style:

   .. code:: matlab

      function c = sum(a, b)
      % Sums two numbers
      %
      % Adds the value of a and b.
      %
      % Returns:
      %   double: The sum of the parameters

   When there are multiple return values, we use another style:

   .. code:: matlab

      function [c d] = sorted2(a, b)
      % Sorts two numbers
      %
      % Returns a and b as c and d so as to always satisfy the condition c <= d.
      %
      % Returns
      % -------
      %   c: double
      %     Smallest number
      %   d: double
      %     Largest number

3. Raises (optional)

   An optional section detailling which errors get raised and under what
   conditions.

   TO BE COMPLETED

4. Warnings (optional)

   An optional section detailling which warnings get raised and under what
   conditions, formatted similarly to Raises.

5. Examples (encouraged)

   A section with explicit commands illustrating as clearly as possible one
   or several ways of calling the function. These commands should be
   written in the doctest format. In partiular, they include the expected
   output.

6. See Also (encouraged)

   An optional section used to refer to related code. This will allow for
   easily browsing the documentation. Fully qualified names should be used
   for class objects and methods.

7. Notes (optional)

   This optional section can provide various additional information about
   the code of interest to the user, such as a discussion about the
   algorithm used by the function. Depreciation warnings can also be
   specified here. We haven't yet specified a formal syntax for those. The
   numpydoc convention uses a Sphinx directive which we tend to avoid.

   Before the section, the comment should be broken by a single empty line
   without ``%``. This should stop matlab from parsing, and thus allow us
   to use further formatting such as LaTeX. The content from now on would
   then only be presented in the Sphinx API. (TODO: check this)

8. References (optional)

   (TODO: should we allow references for the **Notes** section here? see
   numpydoc guide)

Documenting functions
---------------------

The documentation comment is given right after the ``function``
declaration.

Documenting classes
-------------------

We document classes immediately after the ``classdef`` declaration.
However, this class document does not address class properties or
constructor parameters.

Class properties are documented by *not* adding a semicolon ``;`` to
each property, and following them by a comment as in below:

.. code:: matlab

    properties (SetAccess = protected)
        group % (`replab.Group`): Group being representation
        field % ({'R', 'C'}): 'R' for a representation on a real vector space, 'C' for a representation on a complex vector space
        dimension % (integer): Representation dimension
    end

Property types are specified before a semicolon, in parenthesis, as for arguments.

Methods are documented as standalone functions. Do not include ``self`` in the list of parameters.
The constructor is documented separately as any method.

Abstract methods should have a single code line in their body
``error('Abstract');``.

Method blocks can be assigned to a *method group*, using the following syntax:

.. code:: matlab

   methods % Group operations

      function z = compose(self, x, y)

      end

   end

The group name ``Implementation`` or ``Implementations`` will be ignored: it concerns a block that implements abstract methods. The substring ``(Abstract)`` will be erased from group names: it can be used to mention that methods in a block need to be implemented in a subclass.

Documenting types
-----------------

- row cell vectors of a particular type: ``cell{1,:} of TYPE`` or ``cell{1,n} of TYPE``

- row double vectors ``double(1,n)`` or ``double(1,:)``

- square matrices ``double(n,n)``

- matrices ``double(m,n)`` if `m`, `n` are used elsewhere, or `double(:,:)` if the sizes are not referenced elsewhere

- ``charstring`` for row char vectors used as strings (as opposed to the new ``string`` Matlab type unsupported by Octave)

- ``integer`` for integers stored as doubles

- RepLAB own types using Sphinx reference syntax with backticks

Other points
------------

-  Equations: we allow the use of LaTeX equations within ``$`` delimitations,
   as the numpydoc guide stipulates (to be used reasonably since equations
   will not be interpreted in the REPL).

-  We use sparingly the reST conventions for italics, bold and
   monospace, but not for variable names. Only package, method, function
   and class names should be written between single backticks, so that
   the reference is valid for our infrastructure. One cannot use single
   backticks for argument names. Matlab/YALMIP objects can be referred
   with single backticks, as long as they are in the global scope.

-  Use fully qualified references in the first line comment, because that first line comment
   will be used in descriptions of subclasses in other packages.
