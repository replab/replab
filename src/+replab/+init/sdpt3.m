classdef sdpt3 < replab.init.ExternalDependency

    methods

        function self = sdpt3
            self@replab.init.ExternalDependency('sdpt3', 'sdpt3.m');
        end

        function res = works(self)
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
            res = exist('sdpt3.m') == 2;
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
                       'compiler is available on the system (e.g. package ''liboctave-dev'' not installed', char(10), ...
                       'when using octave). Note that other SDP solvers listed in ', char(10), ...
                       'https://yalmip.github.io/allsolvers/ can also be used instead of SDPT3. If you', char(10), ...
                       'do not wish to use Semi-definite programming functionalities you can comment out', char(10), ...
                       'SDP dependencies in replab_config.m.']);
            end
        end

    end

end
