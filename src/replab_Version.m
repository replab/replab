classdef replab_Version
% Describes and manipulates RepLAB version identifiers

    properties
        major % integer: Major version index
        minor % integer: Minor version index
        patch % integer: Patch version index
        snapshot % logical: Is the version a snapshot
    end

    methods

        function self = replab_Version(major, minor, patch, snapshot)
            self.major = major;
            self.minor = minor;
            self.patch = patch;
            self.snapshot = snapshot;
        end

        function newVersion = incrementedMajor(self)
            newVersion = replab_Version(self.major + 1, 0, 0, self.snapshot);
        end

        function newVersion = incrementedMinor(self)
            newVersion = replab_Version(self.major, self.minor + 1, 0, self.snapshot);
        end

        function newVersion = incrementedPatch(self)
            newVersion = replab_Version(self.major, self.minor, self.patch + 1, self.snapshot);
        end

        function newVersion = asRelease(self)
            newVersion = replab_Version(self.major, self.minor, self.patch, false);
        end

        function newVersion = asSnapshot(self)
            newVersion = replab_Version(self.major, self.minor, self.patch, true);
        end

        function str = tag(self)
            assert(~self.snapshot, 'Cannot tag snapshots');
            str = sprintf('v%d.%d.%d', self.major, self.minor, self.patch);
        end

        function str = toText(self)
            str = sprintf('%d.%d.%d', self.major, self.minor, self.patch);
            if self.snapshot
                str = [str '-SNAP'];
            end
        end

        function updateVersionFile(self)
            [path, ~, ~] = fileparts([mfilename('fullpath') '.m']);
            filename = fullfile(path, '..', 'replab_version.txt');
            fid = fopen(filename, 'w');
            assert(fid ~= -1, 'Cannot open replab_version.txt');
            fprintf(fid, '%s\n', self.toText);
            status = fclose(fid);
            assert(status == 0, 'Cannot close replab_version.txt properly.');
        end

        function disp(self)
            if self.snapshot
                snap = '-SNAP';
            else
                snap = ' (release)';
            end
            fprintf('Version %d.%d.%d%s\n', self.major, self.minor, self.patch, snap);
        end

        function newVersion = prompt(self, label)
            newVersion = [];
            while isempty(newVersion)
                userInput = input(sprintf('%s [%s]:\n', label, self.toText), 's');
                if isempty(userInput)
                    newVersion = self;
                else
                    newVersion = replab_Version.fromText(userInput);
                end
            end
        end

    end

    methods (Static)

        function v = fromText(txt)
        % Parses a version string and creates a `replab_Version` object
        %
        % In case of a parse error, returns ``[]``
            v = [];
            tokens = regexp(txt, '^\s*(\d)+[.](\d)+[.](\d)+((-SNAP)?)\s*$', 'tokens', 'once');
            if ~isa(tokens, 'cell')
                return
            end
            switch length(tokens)
              case 3
                tokens{4} = '';
              case 4
                % ok!
              otherwise
                % errored
                return
            end
            snapshot = ~isempty(tokens{4});
            major = str2num(tokens{1});
            minor = str2num(tokens{2});
            patch = str2num(tokens{3});
            v = replab_Version(major, minor, patch, snapshot);
        end

        function v = current
            v = replab_Version.fromVersionFile;
        end

        function v = fromVersionFile
        % Create a `replab_Version` object from the ``replab_version.txt`` text file
            [path, ~, ~] = fileparts([mfilename('fullpath') '.m']);
            filename = fullfile(path, '..', 'replab_version.txt');
            txt = fileread(filename);
            txt = strrep(txt, sprintf('\n'), '');
            txt = strrep(txt, sprintf('\r'), '');
            txt = strrep(txt, ' ', '');
            v = replab_Version.fromText(txt);
        end

        function bumpMajor
            replab_Version.fromVersionFile.incrementedMajor.updateVersionFile;
            fprintf('Now on version ', replab_Version.fromVersionFile.toText);
        end

        function bumpMinor
            replab_Version.fromVersionFile.incrementedMinor.updateVersionFile;
            fprintf('Now on version ', replab_Version.fromVersionFile.toText);
        end

        function bumpPatch
            replab_Version.fromVersionFile.incrementedPatch.updateVersionFile;
            fprintf('Now on version ', replab_Version.fromVersionFile.toText);
        end

        function setTo(versionString)
            v = replab_Version.fromText(versionString);
            assert(~isempty(v), 'Invalid version string');
            v.updateVersionFile;
        end

    end

end
