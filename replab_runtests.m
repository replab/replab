function result = replab_runtests(withCoverage)
    % result = replab_runtests([withCoverage])
    %
    % replab_runtests tests the library functionalities
    %
    % When the option 'withCoverage' is set to 1, code coverage data is
    % generated.
    
    if nargin < 1
        withCoverage = 0;
    end
    
    % Make sure we are in the current path
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
    
    % Check the presence of the MOcov library if needed
    if withCoverage == 1
        MOcovInPath = false;
        try
            mocov_get_absolute_path('.');
            MOcovInPath = true;
        catch
        end
        if ~MOcovInPath
            error('The MOcov library was not found. Did you run replab_addpaths?')
        end
    end
    
    % Check the presence of a SDP solver
    SDPSolverInPath = false;
    try
        x = sdpvar(2);
        sol = optimize([x >= 0, trace(x) == 1], 0, sdpsettings('verbose',0));
        SDPSolverInPath = (sol.problem >= 0);
    catch
    end
    if ~SDPSolverInPath
        warning('No SDP working SDP solver found, some tests will fail.');
    end

    % calls the relevant test suite
    if withCoverage == 1
        result = moxunit_runtests('tests/codeCoverageHelperFunction.m', '-verbose', ...
            '-with_coverage', '-cover', '.', '-cover_exclude', 'external', '-cover_exclude', 'tests', '-cover_exclude', 'docs_src', ...
            '-cover_json_file', 'coverage.json', '-cover_xml_file', 'coverage.xml', '-cover_html_dir', 'coverage_html');
    else
        result = moxunit_runtests('tests', '-verbose', '-recursive');
    end
    
    % return to the previous path
    cd(initialPath);
end
