% # Natural representation
%
% This document shows how to construct the natural representation of a permutation group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group definition
%
% Consider a permutation group acting on $n=5$ elements

S5 = replab.S(5);
group = S5.subgroup({[2 1 4 5 3]})

% The natural representation is 

rep = group.naturalRep

% The representation of group elements are the expected $5 \times 5$ matrices:

rep.image([2 1 5 3 4])

