function adj = edge2adj(edges, n, weights)
% adjacency matrix of the graph
%
% Returns the adjacency matrix for a graph described by a list
% of edges.
%
% Args:
%     edges (integer (\*,2)): array of vertices linked by an edge
%     n (integer): number of vertices
%     weights (double (\*,1), optional): weight associated to each edge
%
% Returns:
%     adj: adjacency matrix of the graph

    if nargin < 3
        weights = 1;
    end

    if isempty(edges)
        edges = zeros(0,2);
    end

    assert(n >= max(max(edges)));
    assert(isequal(size(edges,2), 2));
    assert(length(size(edges)) <= 2);
    assert(isequal(size(weights,2), 1));
    assert((size(weights,1) == size(edges,1)) || (size(weights,1) == 1));

    edges2way = unique([edges; edges(:, [2 1])], 'rows');
    adj = sparse(edges2way(:,1), edges2way(:,2), weights, n, n);
end
