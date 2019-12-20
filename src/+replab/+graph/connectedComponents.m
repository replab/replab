function [subsets componentIndex start next] = connectedComponents(edges)
% Identifies connected components of a graph
%
% Performs the burning algorithm on the network described by the
% edges given in pairs.
%
% Note: This function treats vertices numbers as names, hence only assumes
% the existence of vertices which are linked with other vertices.
%
% Args:
%     edges (double array) : nx2 array of vertices linked by an edge
%
% Returns:
% --------
%     subsets: cell array
%         connex components
%     componentIndex: sparse horizontal vector
%         the index of the component to which each vertex belongs
%     start: full horizontal vector
%         identifies a first element within each each component
%     next: sparse horizontal vector
%         returns for each vertex a next vertex belonging to the same
%         connected component. Value is 0 for the last element of the set.
%         Useful to iterateover all elements of a connected component.
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
        start = [];
        next = [];
        return;
    end

    if size(edges,2) ~= 2
        error('List of edges has wrong size');
    end
    
    % We call the best available burning algorithm
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

        componentIndex = sparse(1,a,b);
    end
    
    if nargout >= 3
        start = zeros(1,length(subsets));
        for i = 1:length(subsets)
            start(i) = subsets{i}(1);
        end
    end
    
    if nargout >= 4
        c = zeros(1,length(a)-length(subsets));
        d = zeros(size(c));
        
        co = 0;
        for i = 1:length(subsets)
            for j = 1:length(subsets{i})-1
                co = co + 1;
                c(co) = subsets{i}(j);
                d(co) = subsets{i}(j+1);
            end
        end
        
        next = sparse(1,c,d,1,max(max(edges)));
    end
end
