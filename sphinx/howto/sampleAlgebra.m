% # Sampling a group representation algebra
%
% This document shows how to sample the algebra of a representation in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Definition of the representation
%
% Consider a group representation

S5 = replab.S(5);
group = S5.subgroup({[2 1 4 5 3]})
rep = group.naturalRep

% ## Sampling the algebra
%
% Elements of the group representation can be sampled with the command

rep.sample

% These matrices describe the action of the group elements on $\mathbb{R}^5$.

