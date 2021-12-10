% # Computing the order of a permutation group
%
% This document shows how to compute the order of a permutation group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Permutation group
%
% Consider the following group

G = replab.PermutationGroup.of([2 3 4 5 1 7 8 9 10 6], [2 1 4 3 6 5 8 7 10 9])

% ## Order computation
%
% The order of the group can be computed by simply

G.order

% Once computed, the order is a property of the group

G

