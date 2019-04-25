function result = replab_runtests
    % result = replab_runtests
    %
    % replab_runtests tests the library functionalities
    
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    
    % Check the presence of the MOxUnit library
    MOxUnitInPath = false;
    try
        moxunit_set_path;
        MOxUnitInPath = true;
    catch
    end
    if ~MOxUnitInPath
        error('The MOxUnit library was not found. Did you run replab_addpaths?')
    end
    
    result = moxunit_runtests('tests','-verbose');
    
    cd(initialPath);
end
