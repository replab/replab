classdef Atlas
% An atlas of finite groups

    methods (Static)

        function readFolder(folderPath)
            if nargin < 1
                folderPath = fullfile(replab.globals.replabPath, 'atlas');
            end
            files = dir(folderPath);
            entries = cell(1, 0);
            for i = 1:length(files)
                file = files(i);
                if ~file.isdir && replab.compat.endsWith(file.name, '.json')
                    filename = fullfile(folderPath, file.name);
                    contents = fileread(filename);
                    try
                        entries{1, end+1} = replab.AtlasEntry.parse(contents);
                    catch
                        err = lasterror;
                        fprintf('Error while reading %s\n', filename);
                        fprintf(strjoin(replab.longStr(err), '\n'));
                    end
                end
            end
            replab.globals.atlasEntries(horzcat(replab.globals.atlasEntries, entries));
        end

        function R = recognizeTrivial(G)
        % Recognizes if the given group is the trivial group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x| x = id>``
            R = [];
            if ~G.isTrivial
                return
            end
            entry = replab.AtlasEntry.trivial;
            R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G));
        end

        function R = recognizeDihedral(G)
        % Recognizes if the given group is the dihedral group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x, a| a^n = x^2 = id, x a x^-1 = a>``
            assert(G.order > 2);
            R = [];
            x = [];
            a = [];
            if G.order > 2^53-1
                return
            end
            if mod(G.order, 2) == 1
                return
            end
            n = double(G.order/2);
            if mod(n, 2) == 1
                Zn = G.derivedSubgroup;
                if ~(Zn.isCyclic && Zn.order == n)
                    return
                end
            else
                G1 = G.derivedSubgroup;
                if ~(G1.isCyclic && G1.order == n/2)
                    return
                end
                T = G.rightCosetsOf(G1).transversal;
                good = false;
                for i = 1:4
                    Zn = G1.closure(T{i});
                    if Zn.isCyclic && Zn.order == n
                        good = true;
                        break
                    end
                end
                if ~good
                    return
                end
            end
            x = G.sample;
            while Zn.contains(x)
                x = G.sample;
            end
            if ~(G.elementOrder(x) == 2 && all(cellfun(@(a) G.isIdentity(G.composeAll({x a x a})), Zn.generators)))
                return
            end
            Rcyclic = replab.Atlas.recognizeCyclic(Zn);
            a = Rcyclic.isomorphism.imageElement(Rcyclic.atlasEntry.group.generator(1));
            entry = replab.AtlasEntry.dihedral(n);
            R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x a}));
        end


        function R = recognizeKlein(G)
        % Recognizes if the given group is the Klein four-group and provides the generators according to the standard presentation
        %
        % The standard presentation is ``<x, a| a^2 = x^2 = id, x a x^-1 = a>``
            R = [];
            if G.order ~= 4
                return
            end
            if G.isCyclic
                return
            end
            x = G.generator(1);
            a = [];
            for i = 2:G.nGenerators
                g = G.generator(i);
                if ~G.eqv(x, g)
                    a = g;
                    break
                end
            end
            assert(~isempty(a));
            entry = replab.AtlasEntry.kleinFourGroup;
            R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x a}));
        end

        function R = recognizeSymmetric(G)
        % Recognizes if the given group is the symmetric group and provides the generators according to the standard presentationx
            R = [];
            [n r] = replab.util.unfactorial(G.order);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses.classes;
            entry = replab.AtlasEntry.symmetric(n);
            for i = 1:length(C)
                S = C{i};
                s = S.representative;
                if G.elementOrder(s) == n
                    for j = 1:length(C)
                        T = C{j};
                        if G.elementOrder(T.representative) == 2
                            U = T.elements;
                            for k = 1:length(U)
                                t = U{k};
                                if entry.group.isMorphismByImages(G, 'images', {s t})
                                    if G.subgroup({s, t}).order == G.order
                                        R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {s t}));
                                        return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end


        function R = recognizeAlternating(G)
        % Recognizes the alternating group and returns the generators corresponding to the standard presentation
            R = [];
            [n r] = replab.util.unfactorial(G.order*2);
            if r ~= 0
                return
            end
            n = double(n);
            C = G.conjugacyClasses.classes;
            entry = replab.AtlasEntry.alternating(n);
            for i = 1:length(C)
                S = C{i};
                s = S.representative;
                if G.elementOrder(s) == n-2
                    for j = 1:length(C)
                        T = C{j};
                        if G.elementOrder(T.representative) == 3
                            U = T.elements;
                            for k = 1:length(U)
                                t = U{k};
                                if entry.group.isMorphismByImages(G, 'images', {s t})
                                    if G.subgroup({s, t}).order == G.order
                                        R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {s t}));
                                        return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        function R = recognizeCyclic(G)
        % Recognizes if the given group is the cyclic group and provides the group generator
        %
        % The standard presentation is ``<x| x^n = id>``
            R = [];
            if G.isTrivial
                return
            end
            if G.order > 2^53-1
                return
            end
            n = double(G.order);
            x = [];
            if ~G.isCyclic
                return
            end
            x = G.sample;
            while G.elementOrder(x) ~= G.order
                x = G.sample;
            end
            entry = replab.AtlasEntry.cyclic(n);
            R = replab.AtlasResult(G, entry, entry.group.isomorphismByImages(G, 'images', {x}));
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
            R = replab.Atlas.recognizeTrivial(group);
            if ~isempty(R)
                return
            end
            R = replab.Atlas.recognizeCyclic(group);
            if ~isempty(R)
                return
            end
            R = replab.Atlas.recognizeKlein(group);
            if ~isempty(R)
                return
            end
            R = replab.Atlas.recognizeDihedral(group);
            if ~isempty(R)
                return
            end
            R = replab.Atlas.recognizeSymmetric(group);
            if ~isempty(R)
                return
            end
            R = replab.Atlas.recognizeAlternating(group);
            if ~isempty(R)
                return
            end
            R = [];
        end

    end

end
