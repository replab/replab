classdef Atlas
% An atlas of finite groups

    methods (Static)

        function computeAndAdd(G)
        % Computes the irreps of the given group and write the corresponding JSON file to the atlas
        %
        % If the group is already recognized by RepLAB, an error is thrown.
        %
        % Note that the `.FiniteGroup.recognize` result is cached, so to be recognized,
        % the given group ``G`` must be constructed again in user code for the new atlas
        % entry to be found.
        %
        % Args:
        %   G (`.FiniteGroup`): Finite group not present in the atlas
            assert(~isempty(replab.globals.gapBinaryPath), 'The GAP binary path must be set in replab_config.m');
            assert(isempty(G.recognize), 'The group is already in the atlas');
            filename = replab.atl.writeAtlasForPermutationGroupUsingGAP(G.permutationGroup, true);
            entry = replab.AtlasEntry.fromFile(filename);
            entries = replab.globals.atlasEntries;
            entries{1,end+1} = entry;
            replab.globals.atlasEntries(entries);
        end

        function read()
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
        %   `+replab.FiniteIsomorphism` or []: A result in case of positive identification; or ``[]`` if unrecognized.
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
            R = replab.atl.Cyclic.recognize(group);
            if ~isempty(R)
                return
            end
            R = replab.atl.Dihedral.recognize(group);
            if ~isempty(R)
                return
            end
            R = replab.atl.Symmetric.recognize(group);
            if ~isempty(R)
                return
            end
            R = replab.atl.Alternating.recognize(group);
            if ~isempty(R)
                return
            end
            R = [];
        end

    end

end
