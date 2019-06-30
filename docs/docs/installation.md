---
layout: docs
title: Installation
position: 1
---

# {{page.title}}

In this chapter we cover the basics of starting to use **RepLAB**.
 
## Cloning the library

To use the library, either:

- Download the latest release on [GitHub](https://www.github.com/replab/replab/releases), but this download does not include the MOxUnit unit test framework.

- Clone the library from GitHub using the following commands:

```
git clone https://www.github.com/replab/replab
cd replab
git submodule init
git submodule update external/MOxUnit
```

This will create a folder **RepLAB** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers and the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit).

## Setting up the path

To use the library, the **RepLAB** folder must be added in Matlab or Octave. This can be done through the `replab_addpaths` scripts.

## Additional resources
**RepLAB** uses the [YALMIP](https://yalmip.github.io) interface to solve convex optimization problems.
In order to use **RepLAB** for semidefinite programming, you will need to [install](https://yalmip.github.io/download/) and [configure](https://yalmip.github.io/tutorial/installation/) it on your system. 
At least one [SDP solver](https://yalmip.github.io/allsolvers/), such as [SeDuMi](https://github.com/SQLP/SeDuMi), is also needed.
The proper functioning of YALMIP can be checked with the command `yalmiptest`.

If you do not wish to install these tools on you own, you can alternatively use the YALMIP and SDPT3 submodules that are included in **RepLAB** for testing purposes. 
For this, update all submodules by issuing the command

```
git submodule update
```

from within the replab folder. 

Within MATLAB/Octave, use then the following commands to finish the installation of the solver:

```
replab_addpaths(1)
cd external/SDPT3/
install_sdpt3
```

On Octave, the package `liboctave-dev` is required for this command to succeed.

## Testing

The proper installation of **RepLAB** can be checked by running the test commands:

```
replab_addpaths
replab_runtests
```

This checks the proper working of the whole package.
