function [subsets, componentIndex] = connectedComponents(edges)
% Identifies connected components of a graph
%
% Performs the burning algorithm on the network described by the
% edges given in pairs.
%
% Note: This function treats vertices numbers as names, hence only assumes
% the existence of vertices which are linked with other vertices.
%
% Args:
%     edges (integer (\n,2)) : array of vertices linked by an edge
%
% Returns:
% --------
%     subsets: cell array
%         connex components
%     componentIndex: sparse integer (1,\n)
%         the index of the component to which each vertex belongs
%
% Example:
%     >>> replab.graph.connectedComponents([1 2; 2 6; 3 4]); % a graph with 5 nodes labelled 1, 2, 3, 4, 6
%
% See also:
%     `replab.Partition.connectedComponents`
%     `replab.graph.burningAlgorithmFast`
%     `replab.graph.burningAlgorithm`


    if isempty(edges)
        % trivial case
        subsets = {};
        componentIndex = [];
        return;
    end

    if size(edges,2) ~= 2
        error('List of edges has wrong size');
    end

    % We call the best available burning algorithm
    if ~replab.dispatch('exists', 'replab.graph.burningAlgorithm')
        replab.dispatch('register', 'replab.graph.burningAlgorithm', 'Fast', 500, ...
                        @(edges) replab.graph.burningAlgorithmFast(edges));
        replab.dispatch('register', 'replab.graph.burningAlgorithm', 'Fallback', 0, ...
                        @(edges) replab.graph.burningAlgorithm(edges));
    end
    subsets = replab.dispatch('call', 'replab.graph.burningAlgorithm', edges);

    % If required, we also compute the next outputs
    if nargout >= 2
        a = [subsets{:}];
        b = zeros(size(a));
        co = 0;
        for i = 1:length(subsets)
            b(co+[1:length(subsets{i})]) = i;
            co = co + length(subsets{i});
        end

        componentIndex = full(sparse(1,a,b));
    end

end
