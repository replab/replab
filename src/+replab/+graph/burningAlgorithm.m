function [componentIndex, subsets] = burningAlgorithm(nbVertices, edges)
% burning algorithm on a graph
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This is the matlab fallback implementation.
%
% Args:
%     nbVertices (integer) : number of vertices in the graph
%     edges: nx2 array of vertices linked by an edge
%
% Returns:
% --------
%     componentIndex: integer (1,\n)
%         the index of the component to which each vertex belongs
%     subsets: cell array
%         connex components
%
% Example:
%     >>> replab.graph.burningAlgorithm(6, [1 2; 2 6; 3 4]); % a graph with 6 nodes labelled 1, 2, 3, 4, 5, 6
%
% See also:
%     replab.UndirectedGraph.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithmFast

    uniquesLeft = unique(edges);
    subsets = {};
    componentIndex = zeros(1, nbVertices);
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
            
%             % We could also burn the links that were already used, but
%             % this way of doing so is super slow...
%             pairs = pairs(setdiff(1:size(pairs,1), union(sel1,sel2)),:);
        end
        subsets{co1} = sort(set);
        uniquesLeft = setdiff(uniquesLeft, set);
        componentIndex(set) = co1;
    end
end
