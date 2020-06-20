function str = whichAll(item)
% Implementation of the which function for Octave compatibility
%
% Returns the same results as ``which(item, '-all')``.
%
% Assumes that ``item`` is implemented as a ``.m`` file.
%
% Args:
%   item (charstring): MATLAB function to look for
%
% Returns:
%   cell(\*,1) of charstring: Paths to the code files
    if replab.compat.isOctave
        str = {};
        paths = strsplit(path, pathsep);
        for i = 1:length(paths)
            candidate = fullfile(paths{i}, [item '.m']);
            if exist(candidate) == 2
                str{end+1,1} = candidate;
            end
        end
    else
        str = which(item, '-all');
    end
end
