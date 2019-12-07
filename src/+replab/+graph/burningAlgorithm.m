function subsets = burningAlgorithm(edges)
% subsets = burningAlgorithm(edges)
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This is the matlab fallback implementation.
%
% Args:
%     edges: nx2 array of vertices linked by an edge
%
% Returns:
%     subsets: cell array with connex components
%
% Example:
%     replab.graph.burningAlgorithm([1 2; 2 6; 3 4]) % a graph with 5 nodes labelled 1, 2, 3, 4, 6
%
% See also:
%     replab.Partition.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithmFast

    uniquesLeft = unique(edges);
    nbVertices = length(uniquesLeft);
    subsets = {};
    co1 = 0;
    while ~isempty(uniquesLeft)
        co1 = co1 + 1;
        set = uniquesLeft(1);
        co2 = 0;
        while co2 < length(set)
            co2 = co2 + 1;
            element = set(co2);
            sel1 = find(edges(:,1) == element);
            sel2 = find(edges(:,2) == element);
            newElements = unique([edges(sel1,2); edges(sel2,1)])';
            newElements = setdiff(newElements, set);
            set = [set, newElements];
        end
        subsets{co1} = sort(set);
        uniquesLeft = setdiff(uniquesLeft, set);
    end
    %subsets{:} % To see the result
end
