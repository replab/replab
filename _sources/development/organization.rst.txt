Code organization and release process
=====================================

Directories
-----------

-  ``docs``: compiled documentation ready to be served by GitHub Pages,
   and is a Git worktree with the branch ``gh-pages`` checked out.

-  ``external``: external libraries, either included as a Git submodule
   (if available on a public repository such as Github), or vendored
   (like ``vpi``). Here is the procedure to update a git submodule names `project`:

    1. Enter the submodule directory with `cd external/project`
    2. Pull the latest version of the project with `git pull origin master`
    3. Go back to the `external` folder with `cd ..` and update the commit number in the file `submodules.ini` so that it matches the new status of the submodule (this is needed for the autoinstall feature).
    4. Commit these changes in the main RepLAB repository: `git add project submodules.ini`, `git commit -m "submodule project updated"`

-  ``sphinx``: source files for the executable documentation, to be
   compiled by ``replab_generate.m``,

-  ``src/+replab``: the main RepLAB package, containing the main
   classes/methods,

-  ``src/+replab/+subpackages``: implementation files for different
   submodules, or experimental stuff,

-  ``tests``: tests written using MOxUnit, augmented with our laws test
   framework.

The release process
-------------------

RepLAB has a ``major.minor.patch`` version number.

The ``develop`` branch contains a snapshot, i.e. a version number that ends in ``-SNAP``.

The ``master`` branch contains a release, i.e. a version number that does not end in ``-SNAP``.

The ``gh-pages`` contains the documentation pages.

See `~root.replab_release` for a description of the release process.
