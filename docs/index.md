---
layout: home
title:  "Home"
section: "Home"
position: 0
---

Current version: **{{site.replabVersion}}**.

[![Join the chat at https://gitter.im/denisrosset/replab](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/denisrosset/replab)
[![Travis CI](https://travis-ci.com/replab/replab.svg?branch=master)](https://travis-ci.com/replab/replab) [![codecov](https://codecov.io/gh/replab/replab/branch/master/graph/badge.svg)](https://codecov.io/gh/replab/replab)

**RepLAB** provides tools to study representations of finite groups and decompose them numerically. It is compatible with both MATLAB and Octave.


## A simple example

The construction of permutation groups and their representations can be achieved in a few lines of code, including the decomposition of their representations into irreducible representations. See the [**full explanation**](docs/publish/Example.html) ...

## How to start using RepLAB

See the [**Installation guide**](docs/installation.html).

## Work in progress

**RepLAB** is a work-in-progress. In particular:

- RepLAB only works in double floating-point precision.
- Most representations are handled by decomposing preimages in words over generators, which limits the order of groups that can be handled.
- While RepLAB has a basic implementation of the BSGS construction, it does not offer much to work with permutation groups.

## Why RepLAB?

Because no open source library exists to decomposes arbitrary permutation/monomial representations into irreducible representations over the reals. RepLAB implements numerical methods that perform this decomposition up to machine precision.

That said, other libraries working on the same problem space include:

- The [GAP System 3 package AREP](https://www.gap-system.org/Gap3/Packages3/arep.html) by Sebastian Egner and Markus Püschel. 
- [NCSOStools](http://ncsostools.fis.unm.si/documentation/awbd) includes an implementation of the Murota-Kanno-Kojima-Kojima-Maehara algorithm to decompose matrix *-algebras.


## Documentation and Support

- Chat it up on [Gitter](https://gitter.im/denisrosset/replab).
- Check the [tutorial](docs/installation.html).

## Contributors

RepLAB and the group theory/linear algebra libraries it depends on were written by [Denis Rosset](https://github.com/denisrosset) and [Jean-Daniel Bancal](https://github.com/jdbancal).

RepLAB references in the `/external` directory the following libraries:

- As a Git submodule, the [MOxUnit](https://github.com/MOxUnit/MOxUnit) test framework by Nikolaas N. oosterhof.

- A copy of the [VPI](https://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) big integer library by John D'Errico.

Feedback and suggestions are always welcome. We ask participants to follow the guidelines of the [Typelevel Code of Conduct](https://typelevel.org/conduct.html).

## License

RepLAB is (C) 2018-2019 Denis Rosset, Jean-Daniel Bancal and other collaborators, and licensed under the [Mozilla Public License 2.0](https://github.com/replab/replab/LICENSE).
