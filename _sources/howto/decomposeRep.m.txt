% # Decomposing a representation
%
% This document shows how to decompose a group representation with *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Representation definition
%
% Let us construct a $4$-dimensional representation

S4 = replab.S(4);
group = S4.subgroup({[2 1 3 4], [1 2 4 3]})
rep = group.naturalRep

% The representation is decomposed with the command

dec = rep.decomposition

% This representation contains 3 irreducible components.

