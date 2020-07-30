function [edges, n, weights] = adj2edge(adj)
% List of graph edges
%
% Returns the list of edges described in an adjacency matrix.
%
% Args:
%     adj (integers (\*,\*)): adjacency matrix of the graph
%
% Returns:
%     edges (integers (\*,2)): array of vertices linked by an edge
%     n (integer): number of nodes
%     weights (double (\*,1)): weight associates to each edge

    assert(length(size(adj)) <= 2);

    n = max(size(adj));
    
    [a, b, weights] = find(adj);
    edges = [a, b];
    
    if norm(weights - weights(1)) == 0
        weights = weights(1);
    end
end
