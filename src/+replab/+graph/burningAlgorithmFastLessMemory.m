function componentIndex = burningAlgorithmFastLessMemory(nbVertices, edges)
% Fast implementation of the burning algorithm using a bit less memory
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This tries to call the fast C++ implementation and
% returns `+replab.+util.Unsuccessful` if it didn't manage to do so.
%
% This implementation is optimized for using less memory, it only returns
% the list of orbits.
%
% Args:
%     nbVertices (integer): Number of vertices
%     edges (integer(n,2)): Array of vertices linked by an edge
%
% Returns
% -------
%     componentIndex: integer (1,\n)
%         the index of the component to which each vertex belongs, or
%         `+replab.+util.Unsuccessful` if unsuccessful
%
% Returns
% -------
% orbits:
%   vector with orbit index for each vertex, or `+replab.+util.Unsuccessful` if unsuccessful
%
% Example:
%     >>> % replab.graph.burningAlgorithmFastLessMemory(6, [1 2; 2 6; 3 4]); % a graph with 6 nodes labelled 1, 2, 3, 4, 5, 6
%
% See also:
%     replab.UndirectedGraph.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithm

    
    persistent compiledInterface;
    if isempty(compiledInterface)
        compiledInterface = replab.dialects.Compiled('cpp', 1);
    end

    componentIndex = compiledInterface.call(1, nbVertices, edges);
end
