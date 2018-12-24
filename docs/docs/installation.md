---
layout: docs
title: Installation
position: 1
---

# {{page.title}}

In this chapter we cover the basics of starting to use RepLAB.

## Cloning the library

To use the library, either:

- Download the latest release on [GitHub](https://www.github.com/replab/replab/releases), but this download does not include the MOxUnit unit test frameowkr.

- Clone the library from GitHub using the following commands:

```
git clone https://www.github.com/replab/replab
cd replab
git submodules init
git submodules update
```

This will create a folder **RepLAB** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers and the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit).

## Setting up the path

To use the library, the **RepLAB** folder must be added in Matlab or Octave. This can be done through the `replab_addpaths` scripts.


## Testing

The proper installation of **RepLAB** can be checked by running the test commands:

```
replab_addpaths
replab_runtests
```

This checks the proper working of the whole package.
