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

    if replab.compat.isOctave
        % warm up Octave before testing
        replab.compat.warmup;
    end

    % Parse the arguments
    args = struct('slowtests', true, 'doctests', true, 'notebooks', true, 'withCoverage', false);
    args = replab.util.populateStruct(args, varargin);

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
        try
            addpath('tests');
            % Set parameters
            matlabTestSuite.slowtests(args.slowtests);
            matlabTestSuite.doctests(args.doctests);
            matlabTestSuite.notebooks(args.notebooks);
            % Launch tests
            result = matlabLaunchTestCovering;
        catch
            result = false;
        end
        rmpath('tests');
    else
        result = moxunit_runtests('tests', '-verbose', '-recursive');
    end

    % Remove the tests folder to the path
    rmpath([pathStr '/../tests']);

    % return to the previous path
    cd(initialPath);
end
