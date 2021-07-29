% # Sampling the commutant of a representation
%
% This document shows how to sample the commutant of a representation in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Definition of the commutant
%
% Consider a group representation

S5 = replab.S(5);
group = S5.subgroup({[2 1 4 5 3]})
rep = group.naturalRep

% The commutant is obtained with the command

comm = rep.commutant

% ## Sampling the commutant
%
% Elements of the commutant can then be sampled

commEl = comm.sample

% These matrices are invariant under transformation by the group representations

groupEl = rep.sample

groupEl*commEl*groupEl'

