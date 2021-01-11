classdef existingSdpSolver < replab.init.Dependency
% Checks whether a good and working SDP solver is already available

    methods

        function self = existingSdpSolver
            self.name = 'existingSdpSolver';
        end

        function res = inPath(self)
            res = true;
        end

        function decentSDPSolverInPath = works(self)
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

                % If a decent solver was found, we make sure it can actually
                % solve an SDP (e.g. the license is valid ;-)
                if decentSDPSolverInPath
                    sol = solvesdp(F, x(1,2), sdpsettings('verbose',0));
                    if isempty(sol) || ~isequal(sol.problem, 0)
                        decentSDPSolverInPath = false;
                        if replab.globals.verboseInit >= 2
                            disp(['The solver ', solver.tag, ' was found, but it produced the following error when called']);
                            disp('to solve and SDP:');
                            disp(['    ', sol.info]);
                        end
                    else
                        replab.init.log_(2, ['Working SDP solver ', solver.tag, ' found in path']);
                    end
                end
            catch
            end
        end

        function init(self, folderName)
        end

    end

end
