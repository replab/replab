classdef sdpt3 < replab.init.ExternalDependency

    methods

        function self = sdpt3
            self@replab.init.ExternalDependency('sdpt3', 'sdpt3.m');
        end

        function res = sdpt3_works(self)
            res = false;
            try
                [blk, Avec, C, b, X0, y0, Z0] = randsdp([2 2], [2 2], 2, 2);
                options.printlevel = 0;
                sdpt3(blk, Avec, C, b, options);
                res = true;
            catch
            end
        end


        function res = inPath(self)
            assert(exist('compileinterfacedata') == 2, 'Needs to have YALMIP in the path');
            x = sdpvar(2);
            F = [x >= 0, trace(x) == 1];
            [interfacedata,recoverdata,solver,diagnostic] = compileinterfacedata(F, [], [], [], sdpsettings, 0, 0);
            res = isempty(diagnostic);
            if ~isempty(strfind(upper(solver.tag), 'LMILAB'))
                % If LMILAB was identified as the best solver to solve the problem, this means that no good solver was found.
                res = false;
            end
        end

        function res = works(self)
            res = false;
            x = sdpvar(2);
            F = [x >= 0, trace(x) == 1];
            [interfacedata,recoverdata,solver,diagnostic] = compileinterfacedata(F, [], [], [], sdpsettings, 0, 0);
            if ~isempty(diagnostic)
                return
            end
            if ~isempty(strfind(upper(solver.tag), 'LMILAB'))
                replab.init.log_(2, 'LMILAB is not an appropriate SDP solver');
                return
            end
            sol = solvesdp(F, x(1, 2), sdpsettings('verbose', 0));
            if isempty(sol) || sol.problem ~= 0
                if replab.globals.verboseInit >= 2
                    disp(['The solver ', solver.tag, ' was found, but it produced the following error when called']);
                    disp('to solve and SDP:');
                    disp(['    ', sol.info]);
                    disp('Trying to use the embedded solver instead.');
                end
                return
            end
            res = true;
        end

        function init(self, path)
            addpath(path);
            replab.init.log_(1, 'Adding embedded SDPT3 solver to the path');

            % Now we run install_sdpt3
            compilationSuccessful = false;
            logSDPT3 = '';
            try
                logSDPT3 = evalc('install_sdpt3;');
                if ~isempty(regexp(logSDPT3, 'Looking for existing binaries\.\.\.incomplete set found\.'))
                    logSDPT3 = evalc('install_sdpt3 -rebuild;');
                end
                compilationSuccessful = isempty(regexp(logSDPT3, 'SDPT3 was not successfully installed.'));
            catch
            end

            if compilationSuccessful
                SDPT3InPath = true;
                if ~isempty(regexp(logSDPT3, 'Looking for existing binaries\.\.\.none found; building\.\.\.'))
                    replab.init.log_(1, 'Compiled SDPT3 binaries');
                end
                replab.init.log_(2, logSDPT3);
            else
                disp(logSDPT3);
                error(['An error occured while trying to set up the SDPT3 solver. This can happen if no', char(10), ...
                       'compiler is available on the system. The functionalities of the library related', char(10), ...
                       'to Semi-definite programming will be disabled. To remedy this, you can install', char(10), ...
                       'an SDP solver listed in https://yalmip.github.io/allsolvers/ .']);
            end
        end

    end

end
