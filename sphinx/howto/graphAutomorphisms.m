% # Graph automorphisms
%
% This document illustrated how *RepLAB* can be used to compute the automorphism group of a graph.
%
% ## Preparation
% As always, before using *RepLAB* commands, initialize the library:

addpath([pwd, '/../..']);
replab_init('verbose', 0);

% ## Graph definition
%
% In *RepLAB*, undirected graphs are represented by the class
% replab.UndirectedGraph. A graph can be initialized from a list of edges 
% or from an adjacency matrix. For instance, the graph with adjacency matrix

M = [0 1 1 1
     1 0 0 1
     1 0 0 1
     1 1 1 0];

% is constructed by calling

graph = replab.UndirectedGraph.fromAdjacencyMatrix(M)

% ## Automorphism computation
%
% The automorphisms of a graph are obtained by simply calling

group = graph.automorphismGroup

% This returns the group under which the graph is invariant. In the present case, this group is a Klein group of order 4, which admits 2 generators.

% ## Bonus
%
% Graphs can have colored vertices as well as weighted edges. These are taken into account during the automorphism computation.

