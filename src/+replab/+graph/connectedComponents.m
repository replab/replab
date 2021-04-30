function [componentIndex, subsets] = connectedComponents(nbVertices, edges)
% Identifies connected components of a graph
%
% Performs the burning algorithm on the network described by the
% edges given in pairs.
%
% Note: This function treats vertices numbers as names, hence only assumes
% the existence of vertices which are linked with other vertices. This
% contrasts with the `+replab.UndirectedGraph` class, which includes
% vertices with no edge.
%
% Args:
%     nbVertices (integer) : number of vertices in the graph
%     edges (integer (\n,2)) : array of vertices linked by an edge
%
% Returns:
% --------
%     componentIndex: integer (1,\n)
%         the index of the component to which each vertex belongs
%     subsets: cell array
%         connex components
%
% Example:
%     >>> replab.graph.connectedComponents(6, [1 2; 2 6; 3 4]); % a graph with 6 nodes labelled 1, 2, 3, 4, 6
%
% See also:
%     `replab.UndirectedGraph.connectedComponents`
%     `replab.graph.burningAlgorithmFast`
%     `replab.graph.burningAlgorithm`


    if isempty(edges)
        % trivial case
        subsets = {};
        componentIndex = zeros(1, nbVertices);
        return;
    end

    if size(edges,2) ~= 2
        error('List of edges has wrong size');
    end

    % We want to do a dispatch with two output arguments, so we do it more
    % simply
    if nargout == 1
    	disp('nargout == 1, calling replab.graph.burningAlgorithmFast');
        componentIndex = replab.graph.burningAlgorithmFast(nbVertices, edges);
        if isa(componentIndex, 'replab.DispatchNext')
	    	disp('nargout == 1, calling replab.graph.burningAlgorithm');
            componentIndex = replab.graph.burningAlgorithm(nbVertices, edges);
        end
    else
        [componentIndex, subsets] = replab.graph.burningAlgorithmFast(nbVertices, edges);
    	disp('nargout ~= 1, calling replab.graph.burningAlgorithmFast');
        if isa(componentIndex, 'replab.DispatchNext')
	    	disp('nargout ~= 1, calling replab.graph.burningAlgorithm');
            [componentIndex, subsets] = replab.graph.burningAlgorithm(nbVertices, edges);
        end
    end

end
