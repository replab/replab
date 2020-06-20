function str = findInstancesInPath(item)
% Finds all implementations of a function/class in the current path
%
% Returns the same results as ``which(item, '-all')``, except it only considers
% the current path if it is explicitly present .
%
% Assumes that ``item`` is implemented as a ``.m`` file.
%
% Note: this is a copy of replab.init.findInstancesInPath
%
% Args:
%   item (charstring): MATLAB function to look for
%
% Returns:
%   cell(\*,1) of charstring: Paths to the code files
    str = {};
    paths = strsplit(path, pathsep);
    for i = 1:length(paths)
        if ~isequal(paths{i}, '.') % Octave has '.' in the path
            candidate = fullfile(paths{i}, [item '.m']);
            if exist(candidate) == 2
                str{end+1,1} = candidate;
            end
        end
    end
end
