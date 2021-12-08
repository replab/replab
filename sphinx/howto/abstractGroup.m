% # Defining an abstract group
%
% This document illustrates how to define an abstract group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Group definition
%
% An abstract group can be defined from its presentation. For instance,
% the cyclic group of order 3 can be defined as

group = replab.AbstractGroup.fromPresentation('< x | x^3 = 1 >')

