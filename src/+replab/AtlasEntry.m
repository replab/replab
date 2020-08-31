classdef AtlasEntry < replab.Obj
% Identifies a user-defined group as a standard group present in an atlas

    properties
        group % (`.AbstractGroup`): Group
        characterTable % (`.CharacterTable`): Group character table
        automorphisms % (`.Automorphisms` or ``[]``): Group automorphisms
    end

    methods

        function self = AtlasEntry(group, characterTable)
            self.group = group;
            self.characterTable = characterTable;
        end

    end

    methods (Static)

        function A = parse(s)
            J = replab.util.parseJSON(s);
            G = replab.AtlasEntry.parseGroup(J);
            C = replab.AtlasEntry.parseCharacterTable(J, G);
            A = replab.AtlasEntry(G, C);
        end

        function G = parseGroup(J)
            generatorNames = J.group.generatorNames;
            permutationGenerators = cellfun(@cell2mat, J.group.permutationGenerators, 'uniform', 0);
            relators = J.group.relators;
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
            name = J.name;
            P = replab.PermutationGroup.of(permutationGenerators{:});
            G = replab.AbstractGroup(generatorNames, P, relators, name);
        end

        function C = parseCharacterTable(J, group)
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

end
