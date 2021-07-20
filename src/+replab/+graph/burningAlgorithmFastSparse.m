function subsets = burningAlgorithmFastSparse(edges)
% Fast implementation of the burning algorithm
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This tries to call the fast C++ implementation and
% returns `+replab.+util.Unsuccessful` if it didn't manage to do so.
%
% This implementation is optimized for extremely sparse graphs, i.e. graphs
% with many vertices, but very few having connections. Therefore, this
% function only returns the sets of connected components.
%
% Args:
%     edges (integer(n,2)): Array of vertices linked by an edge
%
% Returns
% -------
% subsets:
%   Cell array with connex components, or `+replab.+util.Unsuccessful` if unsuccessful
%
% Example:
%     >>> % replab.graph.burningAlgorithmFast_sparse([1 2; 2 6; 3 4]); % a graph with 5 connected nodes labelled 1, 2, 3, 4, 6
%
% See also:
%     replab.UndirectedGraph.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithm

    
    persistent compiledInterface;
    if isempty(compiledInterface)
        compiledInterface = replab.dialects.Compiled('cpp', 1);
    end

    subsets = compiledInterface.call(1, edges);
end
