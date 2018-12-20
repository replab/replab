---
layout: docs
title: Installation
position: 1
---

# {{page.title}}

In this chapter we cover the basics of starting to use RepLab.


## Cloning the library

RepLab is available for download on [github](https://www.github.com/denisrosset/replab/releases), but it relies on additional software. The best way to install it is thus through the following git commands

```
git clone https://www.github.com/denisrosset/replab
cd replab
git submodules init
git submodules update
```

This will create a folder **RepLab** with all the necessary code, including the [VPI](https://ch.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) library for large integers and the testing suite [MOxUnit](https://github.com/MOxUnit/MOxUnit).


## Setting up the path

To use the library, the **RepLab** folder must be added in Matlab or Octave. This can be done through the `addpath` command.


## Testing

The proper installation of the can be checked by running the test command

```
replab_runtests
```

This checks the proper working of the whole package.
