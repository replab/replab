% # Group intersection
%
% This document illustrated how *RepLAB* can be used to compute the intersection between two permutation groups
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group definition
%
% In *RepLAB*, a permutation group is defined from its generators. For instance,
% we can define the group which permutes arbitrarily the first three elements among four element with the following two generators:

generators1 = {[2 1 3 4], [2 3 1 4]};
S4 = replab.S(4);
G1 = S4.subgroup(generators1)

% Similarly, the group which permutes any of the (hence leaves 1 invariant)
generators2 = {[1 3 2 4], [1 3 4 2]};
G2 = S4.subgroup(generators2)

% ## Intersection
%
% The intersection between these two groups is then simply obtained by calling

G12 = intersection(G1, G2)

% In this case, the intersection is the group which only permutes the two middle elements, i.e. S(2) action on elements 2 and 3.

