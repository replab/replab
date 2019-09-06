#Critical Octave bugs

Here is a list of Octave bugs that are problematic to the project. Their resolution will allow to improve the library.

##Bugs affecting replab.CommutantVar:

- to enable automatic copy of handle classes: https://savannah.gnu.org/bugs/?51317
- to split the class methods into individual files : https://savannah.gnu.org/bugs/?54941
- to set up the proper class precedence with respect to yalmip objects (due to this bug, addition of replab.CommutantVar and sdpvar is not commutative (!)) : https://savannah.gnu.org/bugs/?56864
- `isscalar(sdpvar(3))` returns `1` on octave: https://savannah.gnu.org/bugs/?44498