% # Symmetric SDPs
%
% This document shows how to define an SDP variable with *RepLAB* which satisfies some symmetry constraints.

% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Invariant SDP matrix
%
% Consider that we wish to describe a $3 \times 3$ matrix that is invariant
% under cyclic relabelling of its rows and columns. In other words, the matrix
% should satisfy the constraint M(permutation, permutation) == M for

permutation = [2 3 1];

% This is achieved by parametrizing the corresponding commutant with the following command

M = replab.CommutantVar.fromPermutations({permutation}, 'symmetric', 'real')

% Here, we additionally require the matrix to be symmetric and real.
% The SDP matrix $M$ involves only 2 variables, corresponding to the only degrees of freedom left by the constraints. Here is their parametrization:

see(M.fullMatrix)

