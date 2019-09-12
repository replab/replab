#Critical Octave bugs

Here is a list of Octave bugs that are problematic to the project. Their resolution will allow to improve the library.

##Bugs affecting replab.CommutantVar:

- [51317](https://savannah.gnu.org/bugs/?51317): to enable automatic copy of handle classes
- [54941](https://savannah.gnu.org/bugs/?54941): to split the class methods into individual files
- [56864](https://savannah.gnu.org/bugs/?56864): to set up the proper class precedence with respect to yalmip objects (due to this bug, addition of replab.CommutantVar and sdpvar is not commutative (!))
- [44498](https://savannah.gnu.org/bugs/?44498): `isscalar(sdpvar(3))` returns `1` on octave

##Bugs affecting other parts

- `isequal` fails in Octave when objects reference themselves; i.e. our `NiceFiniteGroup` base instance has `self.parent.parent == self.parent` and that throws Octave's `isequal` in a loop
