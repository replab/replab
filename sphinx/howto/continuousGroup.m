% # Defining a continuous group
%
% This document illustrates how to define a compact group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group definition
%
% The group of unitary matrices in dimension 2 can be constructed simply as

group = replab.U(2)

% *RepLAB* integrates several standard compact groups, which can be constructed
% in a similar fashion. This includes in particular the special unitary group,
% the (special) orthogonal group and the symplectic group.

