function edges = adj2edge(adj)
% liste of graph edges
%
% Returns the list of (undirected) edges for a graph described by
% a given adjacency matrix.
%
% Args:
%     adj: adjacency matrix of the graph
%
% Returns:
%     edges: nx2 array of vertices linked by an edge
    
    if ~isequal(size(adj,1), size(adj,2)) || (length(size(adj)) > 2)
        error('The adjacency matrix should be of size nxn');
    end
    
    [a, b] = find(adj);
    edges = [a b];
    edges = unique(sort(edges,2), 'rows');
end
