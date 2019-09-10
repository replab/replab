---
layout: docs
title: API Documentation style guide
---

# API Documentation style guide

This document describes the syntax and best practices for the documentation comments embedded in the source code. 

## Overview

This document is heavily inspired by the [numpydoc docstring guide](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard) and follows its structure. 

We employ the variants proposed by the [Google style](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html) for compactness.

We aim to achieve the following middle ground:

- The documentation comments should be readable on the MATLAB command line when the `help FUN` command is used.

- The documentation comments can be parsed with Sphinx to produce stylish and naviguable documentation.

## Documentation comments

MATLAB documentation call the documentation strings that appears at the beginning of a class/function/method/property definition 'comments'. To distinguish them from code comments, we call those special comments 'documentation comments`, noting that the `docstring` terminology is specific to the Python language (standalone strings in the source code).

The documentation comments are written using a subset of the [reStructuredText](http://docutils.sourceforge.net/rst.html) syntax. This style guide does not list which syntax elements are permitted: use common sense to employ the conventions that are compatible with reStructuredText, and display well standalone on the terminal (as is advocated by the numpydoc guide linked above).

In particular, we use single backticks for class, function, method, property names, and double backticks for verbatim text.

> A guiding principle is that human readers of the text are given precedence over contorting docstrings so our tools produce nice output.

The length of comment lines should be kept to 75 characters; note that this applies to the comment string itself, not including the `%` character, whitespace, and source code tokens. RepLAB uses longer source code lines.

## Sections

The documentation comment consists of a number of sections separated by headings. Each heading is given by a keyword from the following list, with the section content indented.

We use the following section headings:

- `Args` as an alias to Sphinx `Parameters`

- `Notes`

- `Returns` (for multiple return arguments, see the use of underlined section headers below).

- `Raises`

- `See Also`

See also the [Napoleon](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) documentation for other possible section headers.

In order, the sections of a funciton comment are:

1. Short summary
   
   A one line summary that does not use variable names of the function name, e.g.
   
   ```matlab
   
   function c = add(a, b)
   % The sum of two numbers
   ```

2. Deprecation warning

   We haven't yet specified a formal syntax for those. The numpydoc convention uses a Sphinx directive which we tend to avoid.
   
3. Extended summary

	A few sentences giving an extended description. This section should be used to clarify *functionality*, not to discuss implementation details or background theory, which should rather be explored in the **Notes** section below. You may refer to the parameters and the function name (using single backticks), but the parameter descriptions belong in the **Parameters** section.
	
4. Arguments/Parameters

   Description of the function arguments, keywords and their respective types, using Google style.
   
   ```matlab
   % Args:
   %   x (type): Description of parameter `x`.
   %   y: Description of parameter `y` (with type not specified)
	```

	If it is not necessary to specify an argument, use `optional`:
	
	```matlab
	%   x (type, optional): optional parameter of type `type`
	```
	
	When a parameter can only assume one of a fixed set of values, those values can be listed in braces. If the argument is optional, the default value appears first:
	
	```matlab
	%   field ({'R', 'C'})
    ```
	
	When two or more input parameters have exactly the same type, shape and description, they can be combined:
	
	```matlab
	%   a, b (double): Elements to sum
	```
	
	Matlab/Octave encodes integer values using double floating point numbers. Our convention is to name those values `integer`. In the rare case a primitive integer type is needed, we write that type precisely (as in `int32` or `uint32`). Big integers have type `vpi`, which is the name of the external library that supports them.
	
	Strings represented as char arrays have type `char` (recent Matlab versions have a new `string` type which should not be confused with).
	
	Function handles have type `function_handle`.

5. Returns

   Explanation of the returned values and their types. We distinguish two cases.
   
   First, when a single value is returned, we use the Google style:
   
   ```matlab
   function c = sum(a, b)
   % Sums two numbers
   %
   % Returns:
   %   double: The sum of the parameters
   ```
   
   
   ```matlab
   function [c d] = sorted2(a, b)
   % Sorts two numbers
   %
   % Returns
   % -------
   %   c: double
   %     Smallest number
   %   d: double
   %     Largest number
   ```
   
6. Raises

	An optional section detailling which errors get raised and under what conditions.
	
	TO BE COMPLETED
	
7. Warns

	An optional section detailling which warnings get raised and under what conditions, formatted similarly to Raises.
	
8. Warnings

	An optional section with cautions to the user in free text/reST.
	
9. See Also

	An optional section used to refer to related code.
	
	TO BE COMPLETED
	
10. Notes

	An optional section that provides additional information about the code, possibly including a discussion of the
	algorithm.
	
	(TODO: should we allow the use of LaTeX equations with `$` delimitations there, as in the numpydoc guide stipulates?)
	(Is there a way to filter out this section when the `help` function is called from the REPL?)
	
11. References

    (TODO: should we allow references for the **Notes** section here? see numpydoc guide)

12. Examples

	(TODO: allow examples that would be executed by the test framework)

## Documenting functions

The documentation comment is given right after the `function` declaration.

## Documenting classes

We document classes immediately after the `classdef` declaration. However, this class document does not address class properties or constructor parameters.

Class properties are documented by *not* adding a semicolon `;` to each property, and following them by a comment as in below:

```matlab
    properties (SetAccess = protected)
        group % Group being representation
        field % 'R' for a representation on a real vector space, 'C' for a representation on a complex vector space
        dimension % Representation dimension
    end
```

Methods are documented as are standalone functions. Do not include `self` in the list of parameters. The constructor is documented separately as any method.

Abstract methods should have a single code line in their body `error('Abstract');`.

## Other points

- Equations: we do not use currently LaTeX formatting. (TODO: should we, in the future?)

- We use sparingly the reST conventions for italics, bold and monospace, but not for variable names. Variable, package, method, function and class names should be written between single back-ticks.
