% # Defining representation
%
% This document shows how to construct the defining representation of a compact group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group definition
%
% Consider the group of unitary matrices in dimension $d=3$

group = replab.U(3);

% ## Representation
%
% The defining representation is given by

rep = group.definingRep

