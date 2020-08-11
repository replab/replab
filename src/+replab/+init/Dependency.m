classdef Dependency < replab.Str
% Describes a RepLAB dependency
    properties (SetAccess = protected)
        name % (charstring): Name of the dependency
        testFilename % (charstring): File to test the existence of
        inPathFun % (charstring): Function handle that returns true if the library is in the path
        worksFun % (function_handle): Function handle that returns true if the library is working (e.g. by solving a toy problem), if it throws = false, no argument
        initFun % (function_handle): Function handle that initializes the library (e.g. adds it to the path), takes a folder name argument
        githubUsername % (charstring): Name of the GitHub user
        githubRepository % (charstring): Name of the GitHub repository
        gitCommit % (charstring): MD5 of the commit
    end

    methods (Static)


        function d = MOcov
            d = replab.init.Dependency;
            d.name = 'MOcov';
            d.inPathFun = @() any(exist('mocov') == 2);
            d.testFilename = 'MOcov/mocov.m';
            d.worksFun = @() ~isempty(mocov_get_absolute_path('.'));
            d.initFun = @(path) addpath(fullfile(path, 'MOcov'));
            d.githubUsername = 'MOcov';
            d.githubRepository = 'MOcov';
            d.gitCommit = '2078bfc5d28a2cd903767ce5f48c81762645ecca';
        end

        function res = MOxUnit_works
            res = false;
            try
                assertEqual(2 + 2, 4);
                res = true;
            catch
            end
        end

        function d = MOxUnit
            d = replab.init.Dependency;
            d.name = 'MOxUnit';
            d.inPathFun = @() any(exist('moxunit_util_platform_is_octave') == 2);
            d.testFilename = 'MOxUnit/moxunit_set_path.m';
            d.worksFun = @() replab.init.Dependency.MOxUnit_works;
            d.initFun = @(path) run(fullfile(path, 'MOxUnit', 'moxunit_set_path.m'));
            d.githubUsername = 'MOxUnit';
            d.githubRepository = 'MOxUnit';
            d.gitCommit = 'da2a9a31334e9daca52052bae1e59a5980fff606';
        end

        function d = SDP
            d = replab.init.Dependency;
            d.name = 'sdpt3';
            d.inPathFun = @() replab.init.Dependency.sdpSolverInPath;
            d.testFilename = 'sdpt3.m';
            d.worksFun = @() replab.init.Dependency.sdpSolverWorks;
            d.initFun = @(path) replab.init.Dependency.initSDPT3(path);
            d.githubUsername = 'sqlp';
            d.githubRepository = 'sdpt3';
            d.gitCommit = '3b880163eba732c396a1a346dac4a350a259249a';
        end

        function d = YALMIP
            d = replab.init.Dependency;
            d.name = 'YALMIP';
            d.inPathFun = @() exist('yalmiptest') == 2;
            d.testFilename = 'yalmiptest.m';
            d.worksFun = @() ~isempty(yalmip('version'));
            d.initFun = @(path) replab.init.Dependency.initYALMIP(path);
            d.githubUsername = 'yalmip';
            d.githubRepository = 'YALMIP';
            d.gitCommit = 'ac10514a048a6074e63665f069d935e57de1bcd3';
        end

        function res = sdpSolverInPath
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

        function res = sdpSolverWorks
            res = false;
            x = sdpvar(2);
            F = [x >= 0, trace(x) == 1];
            [interfacedata,recoverdata,solver,diagnostic] = compileinterfacedata(F, [], [], [], sdpsettings, 0, 0);
            if ~isempty(diagnostic)
                return
            end
            if ~isempty(strfind(upper(solver.tag), 'LMILAB'))
                replab.init.log(2, 'LMILAB is not an appropriate SDP solver');
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

        function res = sdpt3InPath
            res = false;
            try
                [blk, Avec, C, b, X0, y0, Z0] = randsdp([2 2], [2 2], 2, 2);
                options.printlevel = 0;
                sdpt3(blk, Avec, C, b, options);
                res = true;
            catch
            end
        end

        function initSDPT3(path)
            addpath(path);
            replab.init.log(1, 'Adding embedded SDPT3 solver to the path');

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
                    replab.init.log(1, 'Compiled SDPT3 binaries');
                end
                replab.init.log(2, logSDPT3);
            else
                disp(logSDPT3);
                error(['An error occured while trying to set up the SDPT3 solver. This can happen if no', char(10), ...
                       'compiler is available on the system. The functionalities of the library related', char(10), ...
                       'to Semi-definite programming will be disabled. To remedy this, you can install', char(10), ...
                       'an SDP solver listed in https://yalmip.github.io/allsolvers/ .']);
            end
        end

        function initYALMIP(path)
            subfolders = {'demos' 'extras' 'modules' 'modules/bilevel' 'modules/global' 'modules/moment' 'modules/parametric' 'modules/robust' 'modules/sos' 'operators' 'solvers'};
            addpath(fullfile(path));
            for i = 1:length(subfolders)
                addpath(fullfile(path, subfolders{i}));
            end
        end

    end

    methods % Generic dependency management

        function res = inPath(self)
            res = false;
            try
                f = self.inPathFun;
                res = f();
            catch
                res = false;
            end
        end

        function res = works(self)
            res = false;
            try
                f = self.worksFun;
                res = f();
            catch
                res = false;
            end
        end

        function init(self, folderName)
            f = self.initFun;
            f(folderName);
        end

    end

    methods % GitHub specific

        function u = zipUrl(self)
        % Returns the URL to the GitHub .zip file that contains the dependency files
            assert(~isempty(self.githubUsername));
            assert(~isempty(self.githubRepository));
            u = sprintf('https://github.com/%s/%s/archive/%s.zip', self.githubUsername, self.githubRepository, self.gitCommit);
        end

        function f = zipFilename(self)
        % Returns the .zip filename formatted according to the GitHub conventions
            assert(~isempty(self.githubRepository));
            f = sprintf('%s_%s.zip', self.githubRepository, self.gitCommit);
        end

        function zipDownloadIn(self, folder)
        % Downloads the .zip files containing the dependency in the given folder
            url = self.zipUrl;
            zipPath = fullfile(folder, self.zipFilename);
            if replab.compat.isOctave
                [f, success] = urlwrite (url, zipPath);
                assert(success, 'Download did not work');
            else
                outfilename = websave(zipPath, url);
            end
        end

    end

    methods % RepLAB specific paths

        function preAutoInstallChecks(self)
        % Performs sanity checks: are the directories we are going to modify actually writeable by the current user?
            replabPath = replab.globals.replabPath;
            externalPath = fullfile(replabPath, 'external');

            if exist(externalPath) ~= 7
                error('The external directory %s is not present, but should be by default. Download a fresh copy of RepLAB', externalPath);
            end

            % Test if we have write permission in the parent folder
            [status, attributes] = fileattrib(fullfile(replabPath, '..'));

            if (~status) || (~(attributes.UserWrite || attributes.GroupWrite || attributes.OtherWrite))
                % We cannot read or modify the parent directory
                error('The replab folder %s appears to be located inside a parent folder which cannot be modified by the user. This is not a good idea. Please move it to a folder with write permission enabled and try again.', replabPath);
            end

            % Test if we have write permission in the external folder
            [status, attributes] = fileattrib(externalPath);
            if (~status) || (~(attributes.UserWrite || attributes.GroupWrite || attributes.OtherWrite))
                % We cannot read or modify the external directory
                error('The replab folder %s doesn''t seem to be writeable by the current user. Please fix this and try again', externalPath);
            end
        end

        function autoInstall(self)
        % Installs the dependency automatically in the external/ folder
            self.preAutoInstallChecks;
            replabPath = replab.globals.replabPath;
            externalPath = fullfile(replabPath, 'external');
            depPath = fullfile(replabPath, 'external', self.name);
            zipPath = fullfile(replabPath, 'external', sprintf('%s-%s', self.name, self.gitCommit));
            z = fullfile(replabPath, 'external', self.zipFilename);
            if exist(z) == 2
                replab.init.log(2, 'File %s exists, not downloading.', z);
            else
                replab.init.log(2, 'File %s does not exists, downloading.', z);
                self.zipDownloadIn(fullfile(replab.globals.replabPath, 'external'));
                assert(exist(z) == 2, 'Download failed, file %s did not download.', z);
            end
            switch exist(depPath)
              case 7
                replab.init.log(2, 'Folder %s exists, deleting', depPath);
                replab.compat.rmdirRec(depPath);
              case 0
                replab.init.log(2, 'Folder %s does not exist', depPath);
              otherwise
                error('File %s should be a directory', depPath);
            end
            unzip(z, externalPath);
            movefile(zipPath, depPath);
        end

        function res = externalFolderPath(self)
        % Returns the full path of the folder that contains the dependency in external/
            res = fullfile(replab.globals.replabPath, 'external', self.name);
        end

        function res = inExternal(self)
        % Returns whether the dependency is present in the external/ folder
            replabPath = replab.globals.replabPath;
            verbose = replab.globals.verboseInit;
            path = fullfile(replabPath, 'external', self.name, self.testFilename);
            res = any(exist(path) == [2 3 6 8]);
        end

        function require(self)
        % Initializes the dependency; downloads the dependency under certain conditions
        %
        % This method first checks if the dependency is already present in the path; if so, it verifies that the
        % dependency functions correctly, and throws an error if there is a malfunction.
        %
        % If the dependency is not in the path, but is present in the external/ folder, it initializes it and
        % verifies it works correctly.
        %
        % If the dependency is neither in the path, neither in the external/ folder, it does as follows:
        %
        % * If the `+replab.+globals.autoInstall` flag is ``false``, it throws an error.
        %
        % * If the `+replab.+globals.autoInstall` flag is ``true``, but a ``.git/`` directory is present, it throws an error
        %   mentioning that the Git submodules need to be initialized.
        %
        % * If the `+replab.+globals.autoInstall` flag is ``true`` and no ``.git/`` directory is present, it first checks
        %   whether a ZIP file is present in the ``external/`` folder. If not, it attempts to download it from GitHub.
        %   Then, having the ZIP file, it unzips it and restarts the initialization.
            replabPath = replab.globals.replabPath;
            verbose = replab.globals.verboseInit;
            if ~self.inPath
                replab.init.log(2, 'Dependency %s not present in current path', self.name);
                if ~self.inExternal
                    replab.init.log(2, 'Dependency %s not present in external/ folder', self.name);
                    if replab.globals.autoInstall
                        replab.init.log(1, 'Installing dependency %s', self.name);
                        self.autoInstall;
                    elseif exist(fullfile(replabPath, '.git')) == 7
                        error('Dependency %s not present. Please run ''git submodule init'' and ''git submodule update'' from the command line. Or download the latest RepLAB release without the Git data using the link: https://github.com/replab/replab/archive/master.zip', self.name);
                    else
                        error('Dependency %s not present. Please run replab_init -autoinstall if you have Internet access, or follow installation instructions.', self.name)
                    end
                end
                assert(self.inExternal, 'Auto installation of dependency %s failed', self.name);
                replab.init.log(1, 'Initializing dependency %s', self.name);
                self.init(self.externalFolderPath);
                assert(self.inPath, 'Initialization of dependency %s failed', self.name);
            else
                replab.init.log(2, 'Dependency %s present in current path', self.name);
            end
            if self.works
                replab.init.log(2, 'Dependency %s is working properly', self.name);
            else
                error('Dependency %s available in path, but not working properly', self.name);
            end
        end

    end

end
