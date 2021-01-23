function subsets = burningAlgorithmFast2(nbVertices, edges)
% Fast implementation of the burning algorithm
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This tries to call the fast C++ implementation and
% returns `+replab.DispatchNext` if it didn't manage to do so.
%
% Args:
%     nbVertices (integer): Number of vertices
%     edges (integer(n,2)): Array of vertices linked by an edge
%
% Returns
% -------
% subsets:
%   Cell array with connex components, or `+replab.DispatchNext` if unsuccessful
%
% Example:
%     >>> % replab.graph.burningAlgorithmFast2(6, [1 2; 2 6; 3 4]); % a graph with 6 nodes labelled 1, 2, 3, 4, 5, 6
%
% See also:
%     replab.UndirectedGraph.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithm

    
    persistent compiledInterface;
    if isempty(compiledInterface)
        compiledInterface = replab.dialects.Compiled('cpp', 2);
    end

    [orbits, subsets] = compiledInterface.call(nbVertices, edges);
end
