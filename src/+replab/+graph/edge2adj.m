function adj = edge2adj(edges)
% adj = edge2adj(edges)
%
% Returns the adjacency matrix for a graph described by a list
% of edges.
%
% Args:
%     edges: nx2 array of vertices linked by an edge
%
% Returns:
%     adj: adjacency matrix of the graph
    
    if ~isequal(size(edges, 2), 2) || (length(size(edges)) > 2)
        error('The vector of edges should be of size nx2');
    end
    
    edges2way = unique([edges; edges(:, [2 1])], 'rows');
    adj = sparse(edges2way(:,1), edges2way(:,2), 1);
end
