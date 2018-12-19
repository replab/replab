function result = replab_runstests
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    
    replab_addpaths;
    result = moxunit_runtests('tests','-verbose');
    
    cd(initialPath);
end
