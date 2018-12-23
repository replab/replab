function result = replab_runtests
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    cd(pathStr)
    
    result = moxunit_runtests('tests','-verbose');
    
    cd(initialPath);
end
