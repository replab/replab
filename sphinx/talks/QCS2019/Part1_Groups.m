% ---
% jupyter:
%   jupytext:
%     formats: ipynb,m:light
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.3.2
%   kernelspec:
%     display_name: Octave
%     language: octave
%     name: octave
% ---

% # 1. Constructing groups with RepLAB
%
% We demonstrate a few ways to construct groups using RepLAB. This notebook is a companion to the RepLAB talk at the [Quantum Causal Structures](http://www.cs.ox.ac.uk/conferences/QCS2019/) workshop.

% We start by adding RepLAB to the path. The output below has been generated using MATLAB, but RepLAB is also compatible with recent versions of Octave (>= 4.2).

run ../../../replab_init.m

% ## Group axioms / laws
% RepLAB has knowledge of the group axioms, and can verify them using random sampling of elements. By default, 20 instances are checked for every law. Inspired by [QuickCheck](https://en.wikipedia.org/wiki/QuickCheck).

S4 = replab.S(4);
replab.GroupLaws(S4).check

% ## Permutation groups
% RepLAB uses row vectors of 1-based integers as group elements, which facilitates integration with existing MATLAB code. When constructing a permutation group `replab.S(n)`, we construct an instance of `replab.FiniteGroup` that knows how to interpret those integer row vectors as permutations.

S4 = replab.S(4); % or replab.Permutations(4)
pi = [3 4 1 2]; rho = [3 2 1 4];
sigma = S4.compose(pi, rho)

pi(rho([1 2 3 4])) == sigma([1 2 3 4])

% ### Permutation subgroups
% We can also construct subgroups of S(n) using a list of generators. Example for conditional probability distributions for two measurements settings with two outcomes.

piInput = [3 4 1 2];
piOutput1 = [2 1 3 4];
piOutput2 = [1 2 4 3];
GAlice = replab.S(4).subgroup({piInput piOutput1 piOutput2})

GAlice.order

GAlice.elements

GAlice1 = replab.S(4).subgroup({piInput piOutput1}); % test removing the last generator
GAlice1.order

% ## Signed permutations
% Signed permutation groups are also supported with similar syntax.

sParties = [3 4 1 2];
sFlip = [-1 -2 -3 -4];
sOther = [2 1 3 -4];
GCHSH = replab.SignedPermutations(4).subgroup({sParties sFlip sOther})

GCHSH.order

GCHSH.elements

% ## Wreath product groups
% We construct the wreath product of S(2) on S(2), representing the symetries of "two measurement settings with two outcomes".

Goutcomes = replab.S(2);
Gsettings = replab.S(2);
Gparty = Gsettings.wreathProduct(Goutcomes)

Gparty.elements
