function replab_addpaths(verbose)
% replab_addpaths([verbose])
%
% Sets up the search path in order to enable all functionalities of the
% RepLAB library.
%
% Args:
%     verbose: controls the display level (optional):
%         0: only produce a display in case of error/warning or for
%             critical cases
%         1: informs of the changes made (default value)
%         2: prints full information
%
% Example:
%     replab_addpaths
    
    %% Parameter checking
    if nargin < 1
        verbose = 1;
    end

    
    %% If everything was already set up previously we exit rapidly
    persistent allGood
    if isempty(allGood)
        allGood = false;
    end
    if allGood
        if verbose >= 2
            disp('Exiting because replab_addpaths was already successfully called earlier.');
        end
        return;
    end    

    
    %% Action -- first adding RepLAB itself
    [pathStr, name, extension] = fileparts(which(mfilename));
    pathStr = strrep(pathStr, '\', '/');
    
    % Check if another instance of RepLAB is already in the path
    currentPathStr = strrep(pwd, '\', '/');
    dirName = currentPathStr(find(currentPathStr=='/',1,'last')+1:end);
    cd ..
    AmIMyself = strrep(which('replab_addpaths'), '\', '/');
    cd(dirName);
    if ~isempty(AmIMyself) && ~isequal(AmIMyself, [pathStr, '/', name, extension])
        error(['Another instance of RepLAB in folder ', fileparts(AmIMyself), ' is already in the path.', char(10), ...
            'Use this one or remove it from the path.']);
    else
        if isempty(AmIMyself)
            addpath(pathStr);
            if verbose >= 1
                disp('Adding RepLAB to the path');
            end
        elseif verbose >= 2
            disp('RepLAB is already in the path');
        end
    end

    % We also add the replab package path
    packagePath = strrep(which('replab_version'), '\', '/');
    packagePath = packagePath(1:find(packagePath=='/',1,'last')-1);
    if ~isempty(packagePath) && ~isequal(packagePath, [pathStr, '/src'])
        error(['Another RepLAB package in folder ', packagePath, ' is already in the path.', char(10), ...
            'Use this one or remove it from the path.']);
    else
        if isempty(packagePath)
            packagePath = [pathStr, '/src'];
            addpath(packagePath);
            if verbose >= 1
                disp('Adding RepLAB package to the path');
            end
        elseif verbose >= 2
            disp('RepLAB package is already in the path');
        end
    end
    
    %% VPI
    
    % Making sure the VPI library is in the path and working
    VPIInPath = false;
    try
        VPIInPath = isequal(num2str(vpi('1234567890098765432100123456789')), '    1234567890098765432100123456789');
    catch
    end
    if ~VPIInPath
        if exist([pathStr '/external/vpi/@vpi/vpi.m']) ~= 2
            error(['The VPI library was not found in the folder ', pathStr, '/external/vpi']);
        else
            addpath([pathStr '/external/vpi']);
            if verbose >= 1
                disp('Adding VPI to the path');
            end
            VPIInPath = true;
        end
    elseif verbose >= 2
        disp('VPI is already in the path');
    end
    
    %% MOxUnit
    
    % Making sure MOxUnit is in the path and working
    MOxUnitInPath = false;
    try
        moxunit_set_path;
        MOxUnitInPath = true;
    catch
    end
    if ~MOxUnitInPath
        if exist([pathStr '/external/MOxUnit/MOxUnit/moxunit_set_path.m']) ~= 2
            warning(['The MOxUnit library was not found in the folder ', pathStr, '/external/MOxUnit', char(10), ...
                'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                'It will not be possible to run the tests.']);
        else
            addpath([pathStr '/external/MOxUnit/MOxUnit']);
            moxunit_set_path;
            if verbose >= 1
                disp('Adding MOxUnit to the path');
            end
            MOxUnitInPath = true;
        end
    elseif verbose >= 2
        disp('MOxUnit is already in the path');
    end
    
    %% YALMIP

    % Making sure YALMIP is in the path and working
    YALMIPInPath = false;
    try
        yalmip('version');
        YALMIPInPath = true;
    catch
    end
    if ~YALMIPInPath
        if exist([pathStr '/external/YALMIP/yalmiptest.m']) ~= 2
            warning(['The YALMIP library was found neither in the path nor in the folder ', pathStr, '/external/YALMIP.', char(10), ...
                'Some functionalities of the library might be disabled.']);
        else
            addpath([pathStr '/external/YALMIP']);
            addpath([pathStr '/external/YALMIP/demos']);
            addpath([pathStr '/external/YALMIP/extras']);
            addpath([pathStr '/external/YALMIP/modules']);
            addpath([pathStr '/external/YALMIP/modules/bilevel']);
            addpath([pathStr '/external/YALMIP/modules/global']);
            addpath([pathStr '/external/YALMIP/modules/moment']);
            addpath([pathStr '/external/YALMIP/modules/parametric']);
            addpath([pathStr '/external/YALMIP/modules/robust']);
            addpath([pathStr '/external/YALMIP/modules/sos']);
            addpath([pathStr '/external/YALMIP/operators']);
            addpath([pathStr '/external/YALMIP/solvers']);
            if verbose >= 1
                disp('Adding embedded YALMIP to the path');
            end
            YALMIPInPath = true;
        end
    elseif verbose >= 2
        disp('YALMIP is already in the path');
    end


    %% SDPT3

    % Making sure a woring SDP solver is in the path and working, otherwise
    % tries to add SDPT3
    decentSDPSolverInPath = false;
    SDPT3InPath = false;
    if YALMIPInPath
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
            try
                [blk, Avec, C, b, X0, y0, Z0] = randsdp([2 2], [2 2], 2, 2);
                options.printlevel = 0;
                sdpt3(blk, Avec, C, b, options);
                SDPT3InPath = true;
            catch
            end
            if ~SDPT3InPath
                if exist([pathStr '/external/SDPT3/sdpt3.m']) ~= 2
                    if verbose >= 2
                        warning(['No suitable SDP solver found. In particular, the SDPT3 library was found', char(10), ...
                            'neither in the path nor in the folder ', pathStr, '/external/SDPT3.', char(10), ...
                            'Some functionalities of the library might be disabled.']);
                    end
                else
                    addpath([pathStr '/external/SDPT3']);
                    % Now we run install_sdpt3
                    cd external/SDPT3;
                    install_sdpt3;
                    cd ../..;
                    SDPT3InPath = true;
                end
            elseif verbose >= 2
                disp('SDPT3 is already in the path');
            end
        elseif verbose >= 2
            disp('An SDP solver is already in the path');
        end
    elseif verbose >= 2
        warning('YALMIP not in path, not checking for an SDP solver')
    end

    %% MOcov
    
    % Making sure MOcov is in the path and working
    MOcovInPath = false;
    try
        mocov_get_absolute_path('.');
        MOcovInPath = true;
    catch
    end
    if ~MOcovInPath
        if exist([pathStr '/external/MOcov/MOcov/mocov.m']) ~= 2
            warning(['The MOcov library was not found in the folder ', pathStr, '/external/MOcov', char(10), ...
                'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                'It will not be possible to run tests with code coverage.']);
        else
            addpath([pathStr '/external/MOcov/MOcov']);
            if verbose >= 1
                disp('Adding MOcov to the path');
            end
            MOcovInPath = true;
        end
    elseif verbose >= 2
        disp('MOcov is already in the path');
    end
    
    %% doctest

    % Making sure MOcov is in the path and working
    doctestInPath = false;
    try
        doctest_colors(stdout);
        doctestInPath = true;
    catch
    end
    if ~doctestInPath
        if exist([pathStr '/external/octave-doctest/inst/doctest.m']) ~= 2
            warning(['The octave-doctest library was not found in the folder ', pathStr, '/external/octave-doctest', char(10), ...
                     'Did you run ''git submodule init'' and ''git submodule update''?', char(10), ...
                     'It will not be possible to run doctests.']);
        else
            addpath([pathStr '/external/octave-doctest/inst']);
            if verbose >= 1
                disp('Adding octave-doctest to the path');
            end
            doctestInPath = true;
        end
    elseif verbose >= 2
        disp('octave-doctest is already in the path');
    end

    
    %% If everything was successful, the next call will be quicker
    if VPIInPath && MOxUnitInPath && YALMIPInPath && (decentSDPSolverInPath || SDPT3InPath) && MOcovInPath && doctestInPath
        allGood = true;
    end
    
end
