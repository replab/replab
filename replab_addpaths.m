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


    %% Let us check the Matlab/Octave version
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if ~isOctave
        platform = 'Matlab';
        minimalVersion = '8.6';
        currentVersion = version;
        currentVersion = currentVersion(1:find(currentVersion==' ',1)-1);
    else
        platform = 'Octave';
        minimalVersion = '4.2.2';
        currentVersion = version;
    end
    if ~isLaterVersion(minimalVersion, currentVersion)
        warning(['Current version of ', platform, ' is ', currentVersion, ' but the minimal supported version is ', minimalVersion, '.']);
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
            
            % If a decent solver was found, we make sure it can actually
            % solve an SDP (e.g. the license is valid ;-)
            if decentSDPSolverInPath
                sol = solvesdp(F, x(1,2), sdpsettings('verbose',0));
                if isempty(sol) || ~isequal(sol, 0)
                    decentSDPSolverInPath = false;
                    if verbose >= 2
                        disp(['The solver ', solver.tag, ' was found, but it produced the following error when called']);
                        disp('to solve and SDP:');
                        disp(['    ', sol.info]);
                        disp('Trying to use the embedded solver instead.');
                    end
                end
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
                    compilationSuccessfull = false;
                    cd external/SDPT3;
                    try
                        install_sdpt3;
                        compilationSuccessfull = true;
                    catch
                    end
                    cd ../..;
                    
                    if compilationSuccessfull
                        SDPT3InPath = true;
                    else
                        warning(['An error occured while trying to set up the SDPT3 solver. This can happen if no', char(10), ...
                                 'compiler is available on the system. The functionalities of the library related', char(10), ...
                                 'to Semi-definite programming will be disabled. To remedy this, you can install', char(10), ...
                                 'an SDP solver listed in https://yalmip.github.io/allsolvers/ .']);
                    end
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

    
    %% If everything was successful, the next call will be quicker
    if VPIInPath && MOxUnitInPath && YALMIPInPath && (decentSDPSolverInPath || SDPT3InPath) && MOcovInPath
        allGood = true;
    end
    
end


function ok = isLaterVersion(minimalVersion, currentVersion)
% isLaterVersion - checks if the current version is old enough
%
% ok = isLaterVersion(minimalVersion, currentVersion)
%
% Args:
%     minimalVersion (string) : threshold version
%     currentVersion (string) : actual version
%
% Returns:
%     ok (boolean) : answer
%
% Example:
%     isLaterVersion('4.2.2', '04.3') % true
%     isLaterVersion('4.2.2', '4.2.2.01') % true

    minimalVersion = [minimalVersion, '.'];
    currentVersion = [currentVersion, '.'];
    
    while (sum(minimalVersion == '.') >= 1) && (sum(currentVersion == '.') >= 1)
        minPoint = find(minimalVersion == '.', 1);
        currPoint = find(currentVersion == '.', 1);
        if str2num(minimalVersion(1:minPoint-1)) < str2num(currentVersion(1:currPoint-1))
            ok = true;
            return;
        elseif str2num(minimalVersion(1:minPoint-1)) > str2num(currentVersion(1:currPoint-1))
            ok = false;
            return;
        end
        minimalVersion = minimalVersion(minPoint+1:end);
        currentVersion = currentVersion(currPoint+1:end);
    end
    
    ok = ~((sum(minimalVersion == '.') >= 1) && (sum(currentVersion == '.') < 1));
end
