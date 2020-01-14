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


## How to start using RepLAB?

Have a look at the [**tutorials**](docs/tutorials/tutorials.html)! They have everything to get you started, from installation instructions to hands-on examples.

The documentation of **RepLAB** is organized along 4 directions, following a recent [discussion on software documentation](https://www.divio.com/blog/documentation/):

- [Tutorial](docs/tutorials/tutorials.html): are short hands-on presentations that give you a taste of the goodness of **RepLAB**
- [How-to guides](docs/howto/howto.html): are concise recipes that show you how to achieve a specific goal
- [Topic guides](docs/topic/guides.html): are understanding-oriented presentations that explain the big picture and the key notions on which this software is built
- [Technical reference](docs/reference/reference.html): contains a complete and accurate description of each object of the library

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

RepLAB and the group theory/linear algebra libraries it depends on were initiated by [Denis Rosset](https://github.com/denisrosset) and [Jean-Daniel Bancal](https://github.com/jdbancal). The project has now more [contributors](https://github.com/replab/replab/graphs/contributors).

RepLAB references in the `/external` directory the following libraries:

- As a Git submodule, the [MOxUnit](https://github.com/MOxUnit/MOxUnit) test framework by Nikolaas N. oosterhof.

- A copy of the [VPI](https://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic) big integer library by John D'Errico.

Feedback and suggestions are always welcome. We ask participants to follow the guidelines of the [Typelevel Code of Conduct](https://typelevel.org/conduct.html).

## License

RepLAB is (C) 2018-2019 Denis Rosset, Jean-Daniel Bancal and other collaborators, and licensed under the [Mozilla Public License 2.0](https://github.com/replab/replab/LICENSE).
