---
layout: docs
title: Installation
position: 1
---

# {{page.title}}

In this chapter we cover the basics of starting to use **RepLAB**.
 
## Downloading the library

### Option 1: download the latest RepLAB release, and manually download companion libraries

First, download the latest release on [GitHub](https://www.github.com/replab/replab/releases), but this download does not include the MOxUnit unit test framework.

Then, download and unpack the following libraries in the `extern/` folders, alternatively anywhere on your MATLAB/Octave path, if you want to use the corresponding features.

- to run tests: [MOxUnit](https://github.com/MOxUnit/MOxUnit)
- to test code coverage: [MOcov](https://github.com/MOcov/MOcov)
- to define convex optimization (SDP) problems and run corresponding tests: [YALMIP](https://github.com/yalmip/YALMIP)
- to solve SDP problems and run corresponding tests: [SDPT3](https://github.com/sqlp/sdpt3)

### Option 2: cloning the library (advanced users)

Clone the library from GitHub using the following command:

```
git clone --recursive https://www.github.com/replab/replab
```

which will download the latest `master` version, and update the Git submodules automatically.

This will create a folder **RepLAB** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers, the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit), and support for semidefinite programming.

## Setting up the path

To use the library, the **RepLAB** folder must be added in Matlab or Octave. This can be done through the `replab_addpaths` scripts.

## Additional resources
**RepLAB** uses the [YALMIP](https://yalmip.github.io) interface to solve convex optimization problems.
In order to use **RepLAB** for semidefinite programming, you will need to [install](https://yalmip.github.io/download/) and [configure](https://yalmip.github.io/tutorial/installation/) it on your system. 
At least one [SDP solver](https://yalmip.github.io/allsolvers/), such as [SeDuMi](https://github.com/SQLP/SeDuMi), is also needed.
The proper functioning of YALMIP can be checked with the command `yalmiptest`.

If you do not wish to install these tools on you own, you can alternatively use the YALMIP and SDPT3 submodules that are included in **RepLAB** for testing purposes. 

Within MATLAB/Octave, use the following commands to finish the installation of the solver:

```
replab_addpaths(2,1)
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
