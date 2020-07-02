function res = initMOcov(verbose)
% Adds the MOcov library to the path if it is not yet present
%
% Args:
%   verbose ({0, 1, 2}): Controls the display level
%
% Returns:
%   logical: True if the library is available

    basePath = replab.globals.replabPath;
    res = false;
    try
        mocov_get_absolute_path('.');
        res = true;
    catch
    end
    if ~res
        if exist([basePath '/external/MOcov/MOcov/mocov.m']) ~= 2
            warning(['The MOcov library was not found in the folder ', basePath, '/external/MOcov', char(10), ...
                'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                'It will not be possible to run tests with code coverage.']);
        else
            addpath([basePath '/external/MOcov/MOcov']);
            if verbose >= 1
                disp('Adding MOcov to the path');
            end
            res = true;
        end
    elseif verbose >= 2
        disp('MOcov is already in the path');
    end
end
