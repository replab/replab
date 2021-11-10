function result = replab_runtests(varargin)
% Tests the library functionalities.
%
% See the ``patternsToExclude`` variable below to exclude files from code coverage.
%
% Keyword Args:
%   slowtests (optional, true): Whether to include tests that can possibly
%                               take a long time to run
%   doctests (optional, true): Whether to include doctests
%   notebooks (optional, true): Whether to run the notebook code
%   withCoverage (optional, false): Whether to activate code coverage
%
% Results:
%     logical: True if all tests passed.

    % Parse the arguments
    args = struct('slowtests', true, 'doctests', true, 'notebooks', true, 'withCoverage', false);
    [args, restArgs] = replab.util.populateStruct(args, varargin);    
    
    % Make sure we are in the current path
    initialPath = pwd;
    [pathStr, name, extension] = fileparts(which(mfilename));
    pathStr = strrep(pathStr, '\', '/');
    cd(pathStr)
    cd ..

    % Add the tests folder to the path
    addpath([pathStr '/../tests']);

    % Set test paramters
    ReplabTestParameters.withCoverage(args.withCoverage);
    ReplabTestParameters.onlyFastTests(~args.slowtests);
    if args.slowtests
        replab.Laws.nRuns(20);
    else
        replab.Laws.nRuns(1);
    end

    % Check the presence of the MOxUnit library
    MOxUnitInPath = false;
    try
        moxunit_set_path;
        MOxUnitInPath = true;
    catch
    end
    if ~MOxUnitInPath
        error('The MOxUnit library was not found. Did you run replab_init?')
    end

    % Check the presence of the MOcov library if needed
    if args.withCoverage
        MOcovInPath = false;
        try
            mocov_get_absolute_path('.');
            MOcovInPath = true;
        catch
        end
        if ~MOcovInPath
            error('The MOcov library was not found. Did you run replab_init?')
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

    if args.doctests
        replab_generate('doctests');
    else
        % We clear the doctest folder
        rp = replab.globals.replabPath;
        testRoot = fullfile(rp, 'tests');
        replab.infra.mkCleanDir(testRoot, 'doctests');
    end

    % Create tests for notebooks
    if args.notebooks
        replab_generate('notebooks');
    else
        % We clear the notebook folder
        rp = replab.globals.replabPath;
        testRoot = fullfile(rp, 'tests');
        replab.infra.mkCleanDir(testRoot, 'notebooks');
    end

    % calls the relevant test suite
    if args.withCoverage
        % Here are the files patterns we don't want to include in the
        % coverage monitoring. These are checked by MOcov individually in
        % each subdirectory.
        patternsToExclude = {'replab_*.m', 'callOriginalHelp.m'};

        % We define the test command
        command = 'moxunit_runtests(''tests/codeCoverageHelperFunction.m'', ''-verbose'', ''-with_coverage'', ''-cover'', ''src'', ''-cover_json_file'', ''coverage.json'', ''-cover_xml_file'', ''coverage.xml'', ''-cover_html_dir'', ''coverage_html''';
        for i = 1:numel(patternsToExclude)
            command = [command, ', ''-cover_exclude'', ''', patternsToExclude{i}, ''''];
        end
        command = [command, ')']

        % and call it
        result = eval(command);
    else
        result = moxunit_runtests('tests', '-verbose', '-recursive');
    end

    % Remove the tests folder to the path
    rmpath([pathStr '/../tests']);

    % return to the previous path
    cd(initialPath);
end
