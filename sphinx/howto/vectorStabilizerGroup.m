% # Vector Stabilizer group
%
% This document illustrated how to define groups which leave a vector invariant in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Vector stabilizer group
%
% Given a vector

v = [3 2 3 1 2 5 4];

% in $\mathbb{R}^7$, the most general group acting on the components which leaves the vector invariant can be obtained by first defining the general group permuting all possible components

S7 = replab.S(7)

% and then asking for its stabilization over the vector

Gv = S7.vectorStabilizer(v)

% This group G1 contains all 4 permutations which leave the vector v invariant:

Gv.elements

% ## Vector stabilization over a specific group
%
% Alternatively, it is possible to identify all the permutations which leave a vector invariant and which belong to a specific group. This is achieved by first defining the group of interest, for instance

G1 = replab.PermutationGroup.of([2 1 3:7], [2 3 1 4:7])

% and then stabilize this group over the vector

G1v = G1.vectorStabilizer(v)

% This group only contains the two permutations of Gv which also belong to G1:

G1v.elements

