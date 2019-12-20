function adj = edge2adj(edges, n)
% adjacency matrix of the graph
%
% Returns the adjacency matrix for a graph described by a list
% of edges.
%
% Args:
%     edges (integer matrix): nx2 array of vertices linked by an edge
%     n (integer): number of vertices
%
% Returns:
%     adj: adjacency matrix of the graph
    
    if isempty(edges)
        edges = zeros(0,2);
    end
    
    if ~isequal(size(edges, 2), 2) || (length(size(edges)) > 2)
        error('The vector of edges should be of size nx2');
    end
    
    edges2way = unique([edges; edges(:, [2 1])], 'rows');
    adj = sparse(edges2way(:,1), edges2way(:,2), 1, n, n);
end
