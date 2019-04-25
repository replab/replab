function replab_addpaths
    % Adds directories that are needed to enable all functionalities of the
    % library to the search path.
    
    [pathStr, name, extension] = fileparts(which(mfilename));
    addpath(pathStr);
    addpath([pathStr '/external/vpi']);
    addpath([pathStr '/external/MOxUnit/MOxUnit']);
    addpath([pathStr '/external/MOxUnit/MOxUnit/util']);
    addpath([pathStr '/external/MOcov/MOcov']);
end
