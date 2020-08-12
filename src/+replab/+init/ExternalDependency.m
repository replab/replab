classdef ExternalDependency < replab.init.Dependency
% Describes a RepLAB dependency provided by an external "vendor"

    properties (SetAccess = protected)
        testFilename % (charstring): File witnessing the presence of the library in a folder
        githubUsername % (charstring): Name of the GitHub user
        githubRepository % (charstring): Name of the GitHub repository
        gitCommit % (charstring): MD5 of the commit
    end

    methods

        function self = ExternalDependency(name, testFilename)
            data = replab.init.ExternalDependency.loadConfig(name);
            self.name = name;
            self.testFilename = testFilename;
            if isequal(data.platform, 'GitHub')
                self.gitCommit = data.commit;
                self.githubUsername = data.username;
                self.githubRepository = data.repository;
            end
        end

    end


    methods (Static)

        function data = loadConfig(dependencyName)
            replabPath = replab.globals.replabPath;
            configFile = fullfile(replabPath, 'external', 'modules.ini');
            txt = fileread(configFile);
            lines = strsplit(txt, {'\r', '\n'});
            lines = cellfun(@strtrim, lines, 'uniform', 0);
            mask = cellfun(@(l) isempty(l), lines);
            lines = lines(~mask); % remove empty lines
            [found ind] = ismember(['[' dependencyName ']'], lines);
            lines = lines(ind+1:end);
            nextSection = find(cellfun(@(l) l(1) == '[', lines), 1);
            if ~isempty(nextSection)
                lines = lines(1:nextSection-1);
            end
            data = struct;
            for i = 1:length(lines)
                kv = strsplit(lines{i}, '=');
                assert(length(kv) == 2, 'Error: key-value pair convention');
                k = strtrim(kv{1});
                v = strtrim(kv{2});
                data.(k) = v;
            end
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
                replab.init.log_(2, 'File %s exists, not downloading.', z);
            else
                replab.init.log_(2, 'File %s does not exists, downloading.', z);
                self.zipDownloadIn(fullfile(replab.globals.replabPath, 'external'));
                assert(exist(z) == 2, 'Download failed, file %s did not download.', z);
            end
            switch exist(depPath)
              case 7
                replab.init.log_(2, 'Folder %s exists, deleting', depPath);
                replab.compat.rmdirRec(depPath);
              case 0
                replab.init.log_(2, 'Folder %s does not exist', depPath);
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
                replab.init.log_(2, 'Dependency %s not present in current path', self.name);
                if ~self.inExternal
                    replab.init.log_(2, 'Dependency %s not present in external/ folder', self.name);
                    if replab.globals.autoInstall
                        replab.init.log_(1, 'Installing dependency %s', self.name);
                        self.autoInstall;
                    elseif exist(fullfile(replabPath, '.git')) == 7
                        error('Dependency %s not present. Please run ''git submodule init'' and ''git submodule update'' from the command line. Or download the latest RepLAB release without the Git data using the link: https://github.com/replab/replab/archive/master.zip', self.name);
                    else
                        error('Dependency %s not present. Please run replab_init -autoinstall if you have Internet access, or follow installation instructions.', self.name)
                    end
                end
                assert(self.inExternal, 'Auto installation of dependency %s failed', self.name);
                replab.init.log_(1, 'Initializing dependency %s', self.name);
                self.init(self.externalFolderPath);
                assert(self.inPath, 'Initialization of dependency %s failed', self.name);
            else
                replab.init.log_(2, 'Dependency %s present in current path', self.name);
            end
            if self.works
                replab.init.log_(2, 'Dependency %s is working properly', self.name);
            else
                error('Dependency %s available in path, but not working properly', self.name);
            end
        end

    end

end
