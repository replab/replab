classdef Atlas
% An atlas of finite groups

    methods (Static)

        function read
            folderPath = fullfile(replab.globals.replabPath, 'atlas');
            files = dir(folderPath);
            entries = cell(1, 0);
            for i = 1:length(files)
                file = files(i);
                if ~file.isdir && replab.compat.endsWith(file.name, '.json')
                    try
                        entries{1, end+1} = replab.AtlasEntry.fromFile(file.name);
                    catch
                        err = lasterror;
                        fprintf('Error while reading %s\n', file.name);
                        fprintf(strjoin(horzcat(err.message, replab.longStr(err.stack).'), '\n'));
                    end
                end
            end
            replab.globals.atlasEntries(entries);
        end

        function R = recognize(group)
        % Attempts to identify the given group
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case of positive identification; or ``[]`` if unrecognized.
            entries = replab.globals.atlasEntries;
            for i = 1:length(entries)
                entry = entries{i};
                if entry.canMatch(group)
                    R = entry.match(group);
                    if ~isempty(R)
                        return
                    end
                end
            end
            % TODO
% $$$             R = replab.Atlas.recognizeTrivial(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
% $$$             R = replab.Atlas.recognizeCyclic(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
% $$$             R = replab.Atlas.recognizeKlein(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
% $$$             R = replab.Atlas.recognizeDihedral(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
% $$$             R = replab.Atlas.recognizeSymmetric(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
% $$$             R = replab.Atlas.recognizeAlternating(group);
% $$$             if ~isempty(R)
% $$$                 return
% $$$             end
            R = [];
        end

    end

end
