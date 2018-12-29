function replab_addpaths
    [pathStr, name, extension] = fileparts(which(mfilename));
    addpath(pathStr);
    addpath([pathStr '/external/vpi']);
    addpath([pathStr '/external/MOxUnit/MOxUnit']);
    addpath([pathStr '/external/MOxUnit/MOxUnit/util']);
end
