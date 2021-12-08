% # Symmetrizing a matrix
%
% This document shows how to perform a Reynolds averaging of a matrix in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Matrix to be symmetrized
%
% Let us consider a random matrix

M = rand(5)

% Our aim is to compute the average $\sum_{g\in G} M(g,g)$ over a group $G$,
% where the group elements $g$ act on the rows and on the columns of our matrix.

% ## Group definition
%
% Let us consider the permutation group over which we wish to perform the average

S5 = replab.S(5);
group = S5.subgroup({[2 1 4 5 3]})

% ## Averaging over the group
%
% The averaging is done by projecting our matrix onto the commutant of the group's natural representation

rep = group.naturalRep;
com = rep.commutant;

MP = com.project(M)

% We can check that the matrix is now invariant over permutation by the group elements

el = group.sample
MP(el,el)

