---
layout: docs
title: Installation
subsection: Tutorials
---

# {{page.title}}

In this chapter we cover the basics of starting to use **RepLAB**.

## Downloading the library

The representation theory code of the library is self-contained. Extended features such as *unit tests*, *code coverage* and *convex optimization* make use of external code. Here are two ways of installing the library with the desired set of features to get started. Choose the one which suits you best.

### Option 1: Download the latest RepLAB release, and use our easy install script.


1. Download the latest release on [GitHub](https://www.github.com/replab/replab/releases), which includes the core code.

2. Launch MATLAB/Octave, run the `replab_easyinstall` script in the root folder which will take care of downloading and installing dependencies, which are:

- to run tests: [MOxUnit](https://github.com/MOxUnit/MOxUnit)
- to test code coverage: [MOcov](https://github.com/MOcov/MOcov)
- to define convex optimization (SDP) problems and run corresponding tests: [YALMIP](https://github.com/yalmip/YALMIP)
- to solve SDP problems and run corresponding tests: [SDPT3](https://github.com/sqlp/sdpt3)

3. Go to the "Initializing the library" section below and follow the usage instructions there.


### Option 2: Clone the library

<details><summary>For advanced users only, **click here to unfold**.</summary>
<p>

Clone the library from GitHub using the following command:

```
git clone --recursive https://www.github.com/replab/replab
```

which will download the latest `master` version, and update the Git submodules automatically.

This creates a folder **RepLAB** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers, the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit), and the tools needed for semidefinite programming.

</p>
</details>

## Initializing the library

To use the library, the **RepLAB** folder must be added in Matlab or Octave. Additional paths are also necessary to enable specific functionalities, as mentioned above, and a few variables mush be initialized. This can be done with
```
replab_init
```

This command checks in particular whether an instance of YALMIP is [available](https://yalmip.github.io/download/) and [configured](https://yalmip.github.io/tutorial/installation/) on your system. If this is not the case, the embedded version of yalmip is used. **RepLAB** uses the [YALMIP](https://yalmip.github.io) interface to solve convex optimization problems. The `replab_init` command also ensures that an [SDP solver](https://yalmip.github.io/allsolvers/) such as [SeDuMi](https://github.com/SQLP/SeDuMi) is properly set up. If this is not the case, it activates the embedded SDPT3 solver. The proper installation of a YALMIP instance can be checked with the command `yalmiptest`.

The command `replab_init` should always be used before running any **RepLAB** command. This command only takes some time to run the first time it is called in a MATLAB/Octave session.

## Testing

The proper installation of **RepLAB** can then be checked by running the test commands:

```
replab_runtests
```

This checks the proper working of the whole package (requires the companion packages for test and convex optimization).
