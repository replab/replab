Octave compatibility
====================

.. module:: +replab


We wish to keep RepLAB fully accessible for free. For this reason, the code
is adapted to run within Octave's specificities. For instance, we avoid
calling long chains of constructors which can cause problem. It is also
important to avoid calling 'clear all' within the tests or in the
documentation since this may clear test data while tests are running on
some versions of Octave.

Note that some limitations may also apply on Octave only. For instance, 
Octave doesn't allow for clickable links to be created in the display, 
and it does not support saving/loading class objects. It is also 

Here is a list of Octave bugs that are problematic to the project. Their
resolution will allow to improve the library.

Bugs affecting replab.CommutantVar:
-----------------------------------

-  `51317 <https://savannah.gnu.org/bugs/?51317>`_: to enable automatic
   copy of handle classes
-  `54941 <https://savannah.gnu.org/bugs/?54941>`_: to split the class
   methods into individual files
-  `56864 <https://savannah.gnu.org/bugs/?56864>`_: to set up the
   proper class precedence with respect to yalmip objects (due to this
   bug, addition of replab.CommutantVar and sdpvar is not commutative
   (!))
-  `44498 <https://savannah.gnu.org/bugs/?44498>`_:
   ``isscalar(sdpvar(3))`` returns ``1`` on octave

Bugs affecting other parts
--------------------------

-  ``isequal`` fails in Octave when objects reference themselves; i.e.
   our `replab.NiceFiniteGroup` base instance has
   ``self.parent.parent == self.parent`` and that throws Octave's
   ``isequal`` in a loop
   

