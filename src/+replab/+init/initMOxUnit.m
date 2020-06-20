function res = initMOxUnit(verbose)
% Adds the MOxUnit library to the path if it is not yet present
%
% Returns:
%   logical: True if the library is available
    basePath = replab.globals.replabPath;
    res = false;
    try
        moxunit_set_path;
        res = true;
    catch
    end
    if ~res
        if exist([basePath '/external/MOxUnit/MOxUnit/moxunit_set_path.m']) ~= 2
            warning(['The MOxUnit library was not found in the folder ', basePath, '/external/MOxUnit', char(10), ...
                     'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                     'It will not be possible to run the tests.']);
        else
            addpath([basePath '/external/MOxUnit/MOxUnit']);
            moxunit_set_path;
            if verbose >= 1
                disp('Adding MOxUnit to the path');
            end
            res = true;
        end
    elseif verbose >= 2
        disp('MOxUnit is already in the path');
    end
end
