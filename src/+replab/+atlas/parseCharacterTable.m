        function C = parseCharacterTable(J, group)
            if ~isfield(J, 'characterTable')
                C = [];
                return
            end
            characters = cellfun(@(m) m.', J.characterTable.characters, 'uniform', 0);
            characters = [characters{:}].';
            characters = replab.cyclotomic(characters);
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
                        images{j} = replab.cyclotomic(img);
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
            C = replab.ComplexCharacterTable(group, 'C', conjugacyClasses, characters, args{:});
        end
