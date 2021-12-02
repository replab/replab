% # Averaging a matrix over a continuous group
%
% This document shows how to compute the average of a matrix $M$ over the representation of a compact group in *RepLAB*.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Matrix to be averaged
%
% Let us consider a random matrix

M = rand(4)

% Our aim is to compute the average $\sum_{g\in G} R_g M R_g^\dag$ over a continuous group $G$,
% where $R_g$ is a representation of the group element $g$.

% ## Group and representation definition
%
% Let us consider the group of special orthogonal matrices with determinant equal to 1 in dimension 2

group = replab.SO(2)

% We consider the 4-dimensional tensor representation of this group $R=SO(2)\times SO(2)$

rep = kron(group.definingRep, group.definingRep)


% ## Averaging over the group
%
% The averaging is done by projecting our matrix onto the commutant of the representation

com = rep.commutant;

MP = com.project(M)

% We can check that the matrix is now invariant over actions of the representation

el = rep.sample;
el*MP*el'

