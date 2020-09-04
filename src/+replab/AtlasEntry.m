classdef AtlasEntry < replab.Obj
% Identifies a user-defined group as a standard group present in an atlas

    properties (SetAccess = protected)
        group % (`.AbstractGroup`): Group
        characterTable % (`.CharacterTable` or ``[]``): Group character table
    end

    methods

        function self = AtlasEntry(group, characterTable)
            self.group = group;
            self.characterTable = characterTable;
        end

        function res = canMatch(self, group)
        % Returns whether this atlas entry can match the given group
            res = false;
            gAE = cellfun(@(d) d.abelianInvariants, group.derivedSeries, 'uniform', 0);
            myAE = cellfun(@(d) d.abelianInvariants, self.group.permutationGroup.derivedSeries, 'uniform', 0);
            if ~isequal(gAE, myAE)
                return
            end
            res = true;
            % TODO: tests based on conjugacy classes
        end

        function m = isomorphism(self, group)
            m = self.group.findIsomorphism(group);
        end

        function laws = laws(self)
            laws = replab.laws.AtlasEntryLaws(self);
        end

    end

    methods (Static, Access = protected)

        function s = permToGap(g)
        % Write a GAP System string representation of a RepLAB permutation
        %
        % Args:
        %   g (permutation): Permutation
        %
        % Returns:
        %   charstring: String representation
            s = sprintf('Inverse(PermList([%s]))', strjoin(arrayfun(@num2str, g, 'uniform', 0), ','));
        end

        function G = parseGroup(J)
            generatorNames = J.group.generatorNames;
            permutationGenerators = cellfun(@cell2mat, J.group.permutationGenerators, 'uniform', 0);
            relators = J.group.relators;
            for i = 1:length(relators)
                if isa(relators{i}, 'cell')
                    relators{i} = replab.fp.Letters.print(cell2mat(relators{i}), generatorNames);
                end
            end
            order = J.group.order;
            switch class(order)
              case 'double'
                assert(order <= 2^53, 'Write orders > 2^53 as strings');
                order = vpi(order);
              case 'char'
                order = vpi(order);
              otherwise
                error('Incorrect order type');
            end
            if isfield(J, 'name')
                name = J.name;
            else
                name = [];
            end
            P = replab.PermutationGroup.of(permutationGenerators{:});
            G = replab.AbstractGroup(generatorNames, P, relators, name);
        end

        function C = parseCharacterTable(J, group)
            if ~isfield(J, 'characterTable')
                C = [];
                return
            end
            classesData = J.characterTable.classes;
            classes = cell(1, length(classesData));
            for i = 1:length(classesData)
                cd = classesData{i};
                if isa(cd, 'char')
                    classes{i} = replab.ConjugacyClass.make(group, cd);
                else
                    classes{i} = replab.ConjugacyClass.make(group, group.niceMorphism.preimageElement(cell2mat(cd)));
                end
            end
            conjugacyClasses = replab.ConjugacyClasses(group, classes);
            characters = cellfun(@(m) m.', J.characterTable.characters, 'uniform', 0);
            characters = [characters{:}].';
            characters = replab.cyclotomic.fromStrings(characters);
            args = cell(1, 0);
            if isfield(J.characterTable, 'irreps')
                n = size(characters, 1);
                irreps = cell(1, n);
                id = J.characterTable.irreps;
                for i = 1:n
                    images = id{i};
                    for j = 1:length(images)
                        img = images{j};
                        img = cellfun(@(m) m.', img, 'uniform', 0);
                        img = [img{:}].';
                        images{j} = replab.cyclotomic.fromStrings(img);
                    end
                    irreps{i} = group.repByImages('C', size(images{1}, 1), 'images', images);
                end
                args = horzcat(args, {'irreps', irreps});
            end
            if isfield(J.characterTable, 'irrepNames')
                args = horzcat(args, {'irrepNames', J.characterTable.irrepNames});
            end
            if isfield(J.characterTable, 'classNames')
                args = horzcat(args, {'classNames', J.characterTable.classNames});
            end
            C = replab.CharacterTable(group, conjugacyClasses, characters, args{:});
        end

    end

    methods (Static) % Entries for common finite groups

        function A = trivial
        % Constructs the atlas entry corresponding to the trivial group
            name = 'Trivial group';
            prmGroup = replab.S(1); % this is a legit permutation representation
            generators = cell(1, 0);
            relators = cell(1, 0);
            % Presentation is empty
            ag = replab.AbstractGroup(generators, prmGroup, relators, name);
            classes = ag.conjugacyClasses;
            characters = replab.cyclotomic.eye(1);
            classNames = {'id'};
            irrepNames = {'id'};
            irreps = {ag.trivialRep('C', 1)};
            ct = replab.CharacterTable(ag, classes, characters, 'irreps', irreps, 'classNames', classNames, 'irrepNames', irrepNames);
            A = replab.AtlasEntry(ag, ct);
        end

        function A = dihedral(n)
        % Constructs the atlas entry corresponding to the dihedral group of order 2*n
            assert(n > 2);
            name = sprintf('Dihedral group of order %d', 2*n);
            % Permutation realization
            X = [n:-1:1];
            A = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X, A);
            % Presentation from the groupprops wiki
            % < x, a | a^n = x^2 = 1, x a x^-1 = a^-1 >
            relators = {['a^' num2str(n)] 'x^2' 'x a x^-1 a'};
            ag = replab.AbstractGroup({'x' 'a'}, prmGroup, relators, name);
            ct = replab.ct.DihedralCharacterTable(n);
            assert(ct.group == ag.permutationGroup);
            A = replab.AtlasEntry(ag, ct.imap(ag.niceMorphism.inverse));
        end

        function A = cyclic(n)
        % Constructs the cyclic group of order n
            assert(n >= 2);
            name = sprintf('Cyclic group C(%d) of order %d', n, n);
            % Permutation realization
            X = [2:n 1];
            prmGroup = replab.PermutationGroup.of(X);
            % standard presentation
            % < x | x^n = 1 >
            ag = replab.AbstractGroup({'x'}, prmGroup, {['x^' num2str(n)]}, name);
            ct = replab.ct.CyclicCharacterTable(n);
            assert(ct.group == ag.permutationGroup);
            A = replab.AtlasEntry(ag, ct.imap(ag.niceMorphism.inverse));
        end

        function A = symmetric(n)

        end

    end

    methods (Static) % Generic construction

        function writeJSONforPermutationGroupUsingGAP(G, filename, irreps)
        % Returns the GAP script that outputs the JSON data corresponding to the given group
        %
        % Args:
        %   G (`.PermutationGroup`): Group to compute the information of
        %   filename (charstring): Path to JSON file to write
        %   irreps (logical): Whether to compute irreducible representations
        %
        % Returns:
        %   charstring: GAP System script that computes the relevant data and outputs it in JSON form
            tfile = tempname();
            fid = fopen(tfile, 'wt');
            fprintf(fid, replab.AtlasEntry.getGAPScript(G, irreps));
            fclose(fid);
            [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
            fid = fopen(filename, 'wt');
            fwrite(fid, result);
            fclose(fid);
        end

        function S = getGAPScript(G, irreps)
        % Returns the GAP script that outputs the JSON data corresponding to the given group
        %
        % Args:
        %   G (`.PermutationGroup`): Group to compute the information of
        %   irreps (logical): Whether to compute irreducible representations
        %
        % Returns:
        %   charstring: GAP System script that computes the relevant data and outputs it in JSON form
            if irreps
                scriptName = 'AtlasEntry_irreps.g';
            else
                scriptName = 'AtlasEntry_noirreps.g';
            end
            firstLine = sprintf('G := Group(%s);;', strjoin(cellfun(@(g) replab.AtlasEntry.permToGap(g), G.generators, 'uniform', 0), ', '));
            rest = fileread(fullfile(replab.globals.replabPath, 'src', '+replab', scriptName));
            S = [firstLine char(10) rest];
        end

        function A = forPermutationGroupUsingGAP(G, irreps)
        % Runs GAP System to compute the character table/representation information about a permutation group
        %
        % Args:
        %   G (`.PermutationGroup`): Group to compute the information of
        %   irreps (logical): Whether to compute irreducible representations
        %
        % Returns:
        %   `.AtlasEntry`: The completed atlas for the given group
            tfile = tempname();
            fid = fopen(tfile, 'wt');
            fprintf(fid, replab.AtlasEntry.getGAPScript(G, irreps));
            fclose(fid);
            [status, result] = system([replab.globals.gapBinaryPath ' -q <' tfile]);
            delete(tfile);
            A = replab.AtlasEntry.parse(result);
        end

        function A = parse(s)
        % Parses a JSON string corresponding to an AtlasEntry
        %
        % Args:
        %   s (charstring): JSON string
        %
        % Returns:
        %   `.AtlasEntry`: The parsed entry
            J = replab.util.parseJSON(s);
            G = replab.AtlasEntry.parseGroup(J);
            C = replab.AtlasEntry.parseCharacterTable(J, G);
            A = replab.AtlasEntry(G, C);
        end

    end

end
