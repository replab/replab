function subsets = burning(pairs)
% subsets = burning(pairs)
%
% Performs the burning algorithm on the network described by the
% edges given in pairs.
%
% Args:
%     pairs: nx2 array of vertices linked by an edge
%
% Returns:
%     subsets: cell array with connex components

    persistent triedFastOption fastOptionWorks
    if isempty(triedFastOption) || isempty(fastOptionWorks)
        triedFastOption = false;
        fastOptionWorks = false;
    end
    
    % If possible, we try to use the optimized code
    if (~triedFastOption) || fastOptionWorks
        firstPartWorks = true;
        if ~triedFastOption
            % We try to use the fast option
            % Make sure we are in the current path
            initialPath = pwd;
            try
                % If the program was never compiled, we try to do so
                [pathStr, name, extension] = fileparts(which('replab.graph.burning'));
                pathStr = strrep(pathStr, '\', '/');
                cd(pathStr)
                if exist(['burning_mex.', mexext], 'file') ~= 2
                    mex('-largeArrayDims','burning_mex.cpp')
                end
            catch
                firstPartWorks = false;
            end
            % return to the previous path
            cd(initialPath);
            triedFastOption = true;
        end
        if firstPartWorks
            % The preparation worked, we try call the optimized method
            fastOptionWorks = true;
            try
                subsets = replab.graph.burning_mex(pairs);
            catch
                fastOptionWorks = false;
            end
        end
        if fastOptionWorks
            return;
        end
    end

    % If we end up here, we must use the matlab implementation of the
    % algorithm
    uniquesLeft = unique(pairs);
    subsets = {};
    co1 = 0;
    while ~isempty(uniquesLeft)
        co1 = co1 + 1;
        set = uniquesLeft(1);
        co2 = 0;
        while co2 < length(set)
            co2 = co2 + 1;
            element = set(co2);
            sel1 = find(pairs(:,1) == element);
            sel2 = find(pairs(:,2) == element);
            newElements = unique([pairs(sel1,2); pairs(sel2,1)])';
            newElements = setdiff(newElements, set);
            set = [set, newElements];
            
%             % We could also burn the links that were already used, but
%             % this way of doing so is super slow...
%             pairs = pairs(setdiff(1:size(pairs,1), union(sel1,sel2)),:);
        end
        subsets{co1} = sort(set);
        uniquesLeft = setdiff(uniquesLeft, set);
    end
    %subsets{:} % To see the result
end
