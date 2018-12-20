function result = replab_runtests
    % Tests the library functionalities
    
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    
    replab_addpaths;
    result = moxunit_runtests('tests','-verbose');
    
    cd(initialPath);
end
