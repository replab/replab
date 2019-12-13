function result = replab_runtests(withCoverage, onlyFastTests)
% result = replab_runtests([withCoverage], [onlyFastTests])
%
% replab_runtests tests the library functionalities
%
% Args:
%     withCoverage (boolean): Enable code coverage (optional, default  value is false)
%     onlyFastTests (boolean): Run only a selection of fast tests
%         (optional, default value if false)
%
% Results:
%     result: test results

    if nargin < 1
        withCoverage = false;
    end
    
    if nargin < 2
        onlyFastTests = false;
    end
    
    % Make sure we are in the current path
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    pathStr = strrep(pathStr, '\', '/');
    cd(pathStr)
    cd ..
    
    % Add the tests folder to the path
    addpath([pathStr '/../tests']);
    
    % Set test paramters
    ReplabTestParameters.withCoverage(withCoverage);
    ReplabTestParameters.onlyFastTests(onlyFastTests);
    if onlyFastTests
        replab.Laws.nRuns(1);
    else
        replab.Laws.nRuns(20);
    end
    
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
    decentSDPSolverInPath = false;
    try
        x = sdpvar(2);
        F = [x >= 0, trace(x) == 1];
        [interfacedata,recoverdata,solver,diagnostic] = compileinterfacedata(F, [], [], [], sdpsettings, 0, 0);
        decentSDPSolverInPath = isempty(diagnostic);
        % If LMILAB was identified as the best solver to solve the
        % problem, this means that no good solver was found.
        if ~isempty(strfind(upper(solver.tag), 'LMILAB'))
            decentSDPSolverInPath = false;
        end
    catch
    end
    if ~decentSDPSolverInPath
        warning('No working SDP solver found, some tests will fail.');
    end
    
    % Create doctests
    if ReplabTestParameters.onlyFastTests == 0
        replab_generatedoctests;
    end
    
    % calls the relevant test suite
    if ReplabTestParameters.withCoverage == 1
        result = moxunit_runtests('tests/codeCoverageHelperFunction.m', '-verbose', ...
            '-with_coverage', '-cover', 'src', '-cover_json_file', 'coverage.json', '-cover_xml_file', 'coverage.xml', '-cover_html_dir', 'coverage_html');
    else
        result = moxunit_runtests('tests', '-verbose', '-recursive');
    end
    
    % Remove the tests folder to the path
    rmpath([pathStr '/../tests']);
    
    % return to the previous path
    cd(initialPath);
end
