function res = initSDP(verbose)
% Makes sure a working SDP solver is in the path and working, otherwise tries to add SDPT3
%
% Returns:
%   logical: True if a decent SDP solver is available
    decentSDPSolverInPath = false;
    SDPT3InPath = false;
    basePath = replab.globals.replabPath;
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
            if isempty(sol) || ~isequal(sol.problem, 0)
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
            if exist([basePath '/external/SDPT3/sdpt3.m']) ~= 2
                if verbose >= 2
                    warning(['No suitable SDP solver found. In particular, the SDPT3 library was found', char(10), ...
                             'neither in the path nor in the folder ', basePath, '/external/SDPT3.', char(10), ...
                             'Some functionalities of the library might be disabled.']);
                end
            else
                addpath([basePath '/external/SDPT3']);
                if verbose >= 1
                    disp('Adding embedded SDPT3 solver to the path');
                end

                % Now we run install_sdpt3
                compilationSuccessfull = false;
                logSDPT3 = '';
                try
                    logSDPT3 = evalc('install_sdpt3;');
                    if ~isempty(regexp(logSDPT3, 'Looking for existing binaries\.\.\.incomplete set found\.'))
                        logSDPT3 = evalc('install_sdpt3 -rebuild;');
                    end
                    compilationSuccessfull = true;
                catch
                end

                if compilationSuccessfull
                    SDPT3InPath = true;
                    if (verbose == 1) && ~isempty(regexp(logSDPT3, 'Looking for existing binaries\.\.\.none found; building\.\.\.'))
                        disp('Compiled SDPT3 binaries');
                    elseif verbose >= 2
                        disp(logSDPT3);
                    end
                else
                    disp(logSDPT3);
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
    res = decentSDPSolverInPath || SDPT3InPath;
end
