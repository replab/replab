Code organization in RepLAB
===========================

Directories
-----------

-  ``docs``: compiled documentation ready to be served by GitHub Pages,

-  ``external``: external libraries, either included as a Git submodule
   (if available on a public repository such as Github), or vendored
   (like ``vpi``),

-  ``sphinx``: source files for the executable documentation, to be
   compiled by ``replab_generate.m``,

-  ``src/+replab``: the main RepLAB package, containing the main
   classes/methods,

-  ``src/+replab/+subpackages``: implementation files for different
   submodules, or experimental stuff,

-  ``tests``: tests written using MOxUnit, augmented with our laws test
   framework.
