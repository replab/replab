function subsets = burningAlgorithmFast(edges)
% subsets = burningAlgorithmFast(edges)
%
% Performs the burning algorithm on the network described by the
% edges given in pairs. This tries to call the fast c++ implementation and
% throws a replab:dispatch:tryNext error if it didn't manage to do so.
%
% Args:
%     edges: nx2 array of vertices linked by an edge
%
% Returns:
%     subsets: cell array with connex components
%
% Example:
%     replab.graph.burningAlgorithmFast([1 2; 2 6; 3 4]) % a graph with 5 nodes labelled 1, 2, 3, 4, 6
%
% See also:
%     replab.Partition.connectedComponents
%     replab.graph.connectedComponents
%     replab.graph.burningAlgorithm

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
                [pathStr, name, extension] = fileparts(which('replab_addpaths'));
                pathStr = strrep(pathStr, '\', '/');
                pathStr = [pathStr, '/src/+replab/+graph'];
                cd(pathStr)
                
                needToCompile = true;
                if (exist(['burningAlgorithmFast_mex.', mexext], 'file') == 2) || (exist(['burningAlgorithmFast_mex.', mexext], 'file') == 3)
                    % If the program was already compiled, we check that the
                    % compilation is up to date
                    if exist('burningAlgorithmFast_timestamp.mat', 'file')
                        savedData = load('burningAlgorithmFast_timestamp.mat');
                        fileProperties = dir('burningAlgorithmFast_mex.cpp');
                        if isequal(savedData.timestamp, fileProperties.date)
                            needToCompile = false;
                        end
                    end
                end
                
                if needToCompile
                    % If the program was never compiled, or the source was
                    % modified, we try to compile it
                    mex('-largeArrayDims','burningAlgorithmFast_mex.cpp')

                    % We save the timestamp corresponding to this
                    % compilation
                    fileProperties = dir('burningAlgorithmFast_mex.cpp');
                    timestamp = fileProperties.date;
                    save('burningAlgorithmFast_timestamp.mat', 'timestamp')
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
                subsets = replab.graph.burningAlgorithmFast_mex(edges);
            catch
                fastOptionWorks = false;
            end
        else
            fastOptionWorks = false;
        end
    end

    if ~fastOptionWorks
        % Inform the dispatcher that this method did not succeed
        error('replab:dispatch:tryNext', 'try next');
    end
end
