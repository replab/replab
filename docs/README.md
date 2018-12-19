QDimSum
=======

This package is maintained by Denis Rosset, Armin Tavakoli and Marc-Olivier Renou. The companion paper is available on the [arXiv](https://arxiv.org/search/quant-ph?searchtype=author&query=Rosset%2C+D) (note: **paper will be on arXiv the week of July 23**).

The project is licensed under the [3-Clause BSD License](https://github.com/denisrosset/qdimsum/blob/master/LICENSE.txt). If you use that software to produce research outcomes, we'd really appreciate a citation to the paper mentioned above.

Installation
------------

1. Verify that you have a recent Matlab version (any version `>= 2015b` should work).

2. Download and install [YALMIP](https://yalmip.github.io/download/).

3. Download and install the SDP solver of your choice. We recommend [MOSEK](https://www.mosek.com/downloads/) (free for academia), or [SeDuMi](http://sedumi.ie.lehigh.edu/) (open source).

4. Run `yalmiptest` and verify that SDP problems can be solved in your configuration.

5. Download the software package `QDimSum` from [this URL](https://github.com/denisrosset/qdimsum/archive/master.zip).

6. Unpack the `qdimsum-master.zip` file in a folder of your choice.

7. Add your chosen folder location (`/qdimsum-master`) to your Matlab path (only that folder and not its subfolders).

8. Run the script [TestRAC22.m](https://github.com/denisrosset/qdimsum/blob/master/TestRAC22.m) present in the `qdimsum-master` folder. The script should finish without any error.

You are now ready to start learning about `QDimSum`.

Simplest example (no symmetrisation)
------------------------------------

We now consider the RAC example of our paper for $n=2$ and $d=2$.

For that, we create a file named [RAC22a.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/RAC22a.m), which is a class that extends the [NVProblem](https://github.com/denisrosset/qdimsum/blob/master/NVProblem.m) class provided by `QDimSum`.

```matlab
classdef RAC22a < NVProblem
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self, rank)
        % the first 4 operators are the states, next 2 a projective measurement
        % the last 2 another projective measurement
            dim = 2;
            X = cell(1, 8);
            for i = 1:4
                X{i} = ...
  qdimsum.Random.pureNormalizedDensityMatrix(dim);
            end
            U = qdimsum.Random.unitary(2);
            X{5} = U*[1 0; 0 0]*U';
            X{6} = U*[0 0; 0 1]*U';
            U = qdimsum.Random.unitary(2);
            X{7} = U*[1 0; 0 0]*U';
            X{8} = U*[0 0; 0 1]*U';
        end
        function K = sampleStateKraus(self)
            K = eye(2); % dimension is 2
        end
        function obj = computeObjective(self, X, K)
            obj = 0;
            for x1 = 1:2
                for x2 = 1:2
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        rho = X{x1+(x2-1)*2};
                        M = X{4+b+(y-1)*2};
                        obj = obj + real(trace(M * rho)/8);
                    end
                end
            end
        end
    end
end
```

Now, a few explanations.

- The property `forceReal = true` indicates to `QDimSum` that we can take the real part of our moment matrix after sampling. All the problems we considered have this property; beware, we have not tested `QDimSum` on problems involving complex moment matrices.

- The method `sampleOperators` is an oracle that returns a generic sample of the operators. In our problem, that involves sampling from the states $\rho_{1,1} \ldots \rho_{2,2}$, and from the measurement operators $M^1_1 \ldots M^2_2$.

- The method `sampleStateKraus` is an oracle that returns a generic sample of the Kraus operator corresponding to a state in a set determined by the problem of interest. For Bell inequalities, this corresponds to the ket of a shared entangled state. For distributed computation tasks, it reduces to the identity operator.

- The method `computeObjective` computes the objective function of the problem of interest for a specific realization given by the operators `X` and the Kraus operator `K`.

Now, to optimize over that problem, we write the following in a [TestRAC22a.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/TestRAC22a.m) script:

```matlab
settings = NVSettings;
problem = RAC22a;
monomials = {'npa' 2};
disp('The objective value is:')
upperBoundSDP = nvOptimize(problem, monomials, ...
  'none', settings)
disp('While the analytical value is:')
sol = 1/2*(1+1/sqrt(2))
```

Again, a few explanations:

- The `settings` variable is an instance of the class [NVSettings](https://github.com/denisrosset/qdimsum/blob/master/NVSettings.m), that parameterizes not only the behaviour of `QDimSum` but also the options passed to YALMIP and the SDP solver. The default options are fine for now.

- The `monomials` variable provides the generating set of monomials used to construct the moment matrix. It can be a cell vector containing list of variable indices, or simply a cell array `{'npa' levelN}` where `levelN >= 1` is an integer corresponding to the maximal degree of monomials in the generating set.

- For now, we pass the parameter `method='none'` as our problem does not define a symmetry group.

- The [nvOptimize](https://github.com/denisrosset/qdimsum/blob/master/nvOptimize.m) function is where the magic happens (see its source code for full documentation about its use).


Using symmetrisation
--------------------

To use symmetrisation, you need to specify the symmetry group of the problem. For that, add the method below to the class `RAC22` (just before or after `computeObjective`). For our purposes, we renamed the class [RAC22b](https://github.com/denisrosset/qdimsum/blob/master/tutorial/RAC22b.m) for that second step.

```matlab
        function generators = symmetryGroupGenerators(self)
            swapX1X2 = [1 3 2 4 7 8 5 6];
            flipX1 = [2 1 4 3 6 5 7 8];
            generators = [swapX1X2
                          flipX1];
        end
```

The `symmetryGroupGenerators` method specifies the symmetries one wishes to exploit. A symmetry corresponds to a permutation of the sequence of physical operators which leaves the objective function invariant. If there are $N$ physical operators in the problem, a symmetry group element is specified by providing the image of the list $(1, \ldots N)$ under that symmetry group element. If the symmetry group has $M$ generators, the method `symmetryGroupGenerators` should return a $M \times N$ matrix containing one such permutation per row.

Here, The generator `swapX1X2` permutes the indices $x_1$ and $x_2$ in $\rho_{x1,x2}$, while `flipX1` permutes $x_1$.

Now, you have access to all the symmetrisation methods of `QDimSum`, we use the following [script](https://github.com/denisrosset/qdimsum/blob/master/tutorial/TestRAC22b.m):

```matlab
settings = NVSettings;
problem = RAC22a;
monomials = {'npa' 2};
nvOptimize(problem, monomials, 'none', settings)
nvOptimize(problem, monomials, 'reynolds', settings)
nvOptimize(problem, monomials, 'isotypic', settings)
nvOptimize(problem, monomials, 'irreps', settings)
nvOptimize(problem, monomials, 'blocks', settings)
```

The application of summarization to a problem amounts to two parts. First we reduce the number of linearly independent basis elements found through the sampling procedure. This number will decrease as the number of applied independent symmetries increases. Second, the final SDP matrix can be block-diagonalized. The size and number of these blocks depends on the symmetry group identified.

The software package offers several methods for symmetrized implementations of the NV hierarchy.

- `none`: No symmetrization; standard implementation of the NV hierarchy.
- `reynolds`: Computes the Reynolds' operator but performs no block-diagonalization of the SDP matrix. 
- `isotypic`: Computes the Reynolds' operator and performs a partial block-diagonalization of the SDP matrix corresponding a decomposition into isotypic components.
- `irreps`: Computes the Reynold’s operator and performs a complete block-diagonalization of the SDP matrix.
- `blocks`: Identifies a complete block-diagonalization and computes the moment matrix directly in the block-diagonal form.

Optionally, if the isotypic or irreducible decomposition of the group action is already known, it can be passed to the `nvOptimize` method (see its [documentation](https://github.com/denisrosset/qdimsum/blob/master/nvOptimize.m)).

Each of these methods offers a different degree of symmetrization which may be more or less relevant depending on the specific problem of interest. In particular, due to a temporary limitation of our code, `irreps` and `blocks` require all the irreducible representations appearing to be of real type.

Note: we allow permutations to flip the sign of the operator. For that, just flip the sign of the integer in the image (see [I3322.m](https://github.com/denisrosset/qdimsum/blob/master/I3322c.m) for an example).

Controlling monomials: families of operators
--------------------------------------------

Instead of the [NPA type hierarchy](http://iopscience.iop.org/article/10.1088/1367-2630/10/7/073013/meta), one may want to use [local levels](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.030501) to decide on the monomials involved.

For that purpose, one has to provide an additional method `operatorTypes`:

```matlab
        function types = operatorTypes(self)
            types = {1:4 5:8};
        end
```

where `types` is a cell array marking each physical operator (states, measurements etc.) that appear in the physical problem of interest. For example, in a prepare and measure scenario with `N` preparations and M measurement operators, the cell array reads `{1:N N+(1:M)}`, i.e., the first `N` positions will correspond to preparations and the final `M` to the measurement operators. Similarly, in a bipartite Bell experiment with `N` measurement operators for each party, the cell array reads `{1:N N+(1:N)}`.

Now, one can specify that the monomials involved are $1$, $\rho$, $M$, $\rho \cdot M$ by running the code as follows:
```matlab
settings = NVSettings;
problem = RAC22c;
monomials = {'families' [] [1] [2] [1 2]};
nvOptimize(problem, monomials, 'blocks', settings)
```

See the files [RAC22c.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/RAC22c.m) and [TestRAC22c.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/TestRAC22c.m).

Controlling settings
--------------------

We describe below the most important options in [NVSettings](https://github.com/denisrosset/qdimsum/blob/master/NVSettings.m). To change them from the default options, one calls `NVSettings` as such:

```matlab
settings = NVSettings( ...
  'yalmipSettings', ...
  sdpsettings('solver', 'sedumi', 'verbose', 0), ...
  'verbosityLevel', 0);
```

This selects SeDuMi as a solver (by passing the relevant options to YALMIP), and prevents the code to output its progress, by silencing both `QDimSum` (with `verbosityLevel`) and SeDuMi (`verbose = 0`).

The documentation of all options can be found in the [NVSettings.m source file](https://github.com/denisrosset/qdimsum/blob/master/NVSettings.m).

Sanity checks
-------------

In the options above, the setting `checkLevel` is set to `1` by default. Then, `QDimSum` will perform consistency checks at every step of the computation. Those checks are not expensive and should not be disabled without good reason. Moreover, one can help `QDimSum` in checking the user provided data by explaining to the toolbox the constraints that the operators have to satisfy.

```matlab
        function C = operatorSDPConstraints(self, X)
            C = X; % all operators are SDP
        end
        function C = operatorEqualityConstraints(self, X)
            dim = 2;
            % X{5}, X{6} form a projective measurement
            % X{7}, X{8} form a projective measurement
            C = {eye(dim) - X{5} - X{6}
                 eye(dim) - X{7} - X{8}};
        end
        function C = scalarEqualityConstraints(self, X)
            dim = 2;
            mm = eye(dim)/dim; % states have the same
			% trace as the maximally mixed state
            C = {mm - X{1}
                 mm - X{2}
                 mm - X{3}
                 mm - X{4}};
        end
```

The documentation of those constraints is given in [NVProblem.m](https://github.com/denisrosset/qdimsum/blob/master/NVProblem.m), here they correspond respectively to: the positive semidefiniteness of the operators, to the projector elements summing to identity, and the states being normalized.

In the file [RAC22d.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/RAC22d.m), we introduced an error on purpose in the symmetry group generators. When we run `TestRAC22d.m`, we get the error message:

```
Error using qdimsum.Check.sampleObeysConstraints (line 68)
For generator 1  3  2  4  7  6  5  8: operator
equality constraint #1 violated, maximal
singular value 0.960222
```

that identifies the problem: run [TestRAC22d.m](https://github.com/denisrosset/qdimsum/blob/master/tutorial/TestRAC22d.m) to reproduce that error and its detection.

Additional tools
----------------

- [findSymmetryGroupGenerators](https://github.com/denisrosset/qdimsum/blob/master/findSymmetryGroupGenerators.m): for a given physical problem and its ambient group (defined using the `ambientGroupGenerators` method, this function finds generators of the corresponding symmetry group. This is a useful tool when the objective function has none or only few *obvious* symmetries found by inspection. Despite the search for symmetries being fairly demanding, it is useful for small scale problems, as well as middle sized problems in which a lesser number of symmetries are already known.

- In [NVSettings](https://github.com/denisrosset/qdimsum/blob/master/NVSettings.m), one controls the number of samples in a batch with `sampleChunkSize`. We also provide the static methods `NVSettings.yalmipMOSEK`, `NVSettings.yalmipSeDuMi`, `NVSettings.yalmipSDPT3`, `NVSettings.yalmipSDPNAL` and `NVSettings.yalmipSCS` that take a single parameter `tolerance` and possible additional parameters. They return a `sdpsettings`-like structure, with the tolerance of the solver set to the provided value, and the possible additional options. To use it as part of the settings, write `settings = NVSettings('verbosityLevel', 0, 'yalmipSettings', NVSettings.yalmipMOSEK(1e-9))`.

- Partial support for rank-constrained problems is available through [RankProblem.m](https://github.com/denisrosset/qdimsum/blob/master/RankProblem.m), but no documentation is currently available.

