% # Group definition from generators
%
% This document illustrates how to define a group from a set of generators in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group from one generator
%
% In *RepLAB*, a generator on $n$ elements is described by a vector in $\mathbb{R}^n$ containing the image of each element. For instance,

genCyclic = [2 3 4 5 1];

% defines a syclic permutation. The group induced by this generator is created by taking the corresponding subgroup of $S(n)$:

S5 = replab.S(5);
groupCyclic = S5.subgroup({genCyclic})

% Note that the generator is inserted in a cell array.

% ## Group from multiple generators
%
% A group generated by several permutations is created similarly as a subgroup of $S(n)$. For instance, the generators

gen1 = [2 1 3 4 5];
gen2 = [2 3 4 5 1];

% generate

groupS5 = S5.subgroup({gen1, gen2})

% where the cell array contains the set of generators.

