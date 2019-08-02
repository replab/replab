---
layout: docs
title: Installation
position: 1
---

# {{page.title}}

In this chapter we cover the basics of starting to use **RepLAB**.
 
## Downloading the library

The representation theory code of the library is self-contained. Extended features such as *unit tests*, *code coverage* and *convex optimization* make use of external code. Here are two ways of installing the library with the desired set of features to get started.

### Option 1: Download the latest RepLAB release, and manually download companion libraries

First, download the latest release on [GitHub](https://www.github.com/replab/replab/releases), which includes the core code.

Then, download and unpack the following libraries in the `extern/` folders if you want to use the corresponding features.

- to run tests: [MOxUnit](https://github.com/MOxUnit/MOxUnit)
- to test code coverage: [MOcov](https://github.com/MOcov/MOcov)
- to define convex optimization (SDP) problems and run corresponding tests: [YALMIP](https://github.com/yalmip/YALMIP)
- to solve SDP problems and run corresponding tests: [SDPT3](https://github.com/sqlp/sdpt3)


### Option 2: Clone the library (advanced users)

Clone the library from GitHub using the following command:

```
git clone --recursive https://www.github.com/replab/replab
```

which will download the latest `master` version, and update the Git submodules automatically.

This creates a folder **RepLAB** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers, the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit), and the tools needed for semidefinite programming.


## Setting up SDP support
**RepLAB** uses the [YALMIP](https://yalmip.github.io) interface to solve convex optimization problems. The steps above make sure the YALMIP and SDPT3 code are downloaded. To activate the SDPT3 solver, run the following command from the `SDPT3` folder
```
install_sdpt3
```
On Octave, the package `liboctave-dev` is required for this command to succeed.


Alternatively, **RepLAB** can also detect automatically any instance of YALMIP which has been [installed](https://yalmip.github.io/download/) and [configured](https://yalmip.github.io/tutorial/installation/) on your system. In particular, any [SDP solver](https://yalmip.github.io/allsolvers/), such as [SeDuMi](https://github.com/SQLP/SeDuMi), can be used by **RepLAB** instead of SDPT3. The proper installation of a YALMIP instance can be checked with the command `yalmiptest`.


## Setting up the path

To use the library, the **RepLAB** folder must be added in Matlab or Octave. This can be done with
```
replab_addpaths
```


## Testing

The proper installation of **RepLAB** can then be checked by running the test commands:

```
replab_runtests
```

This checks the proper working of the whole package (requires the companion packages for test and convex optimization).
