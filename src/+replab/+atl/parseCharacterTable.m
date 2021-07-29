function C = parseCharacterTable(group, field, data)
% Reads a character table from parsed JSON data
%
% Args:
%   group (`+replab.FiniteGroup`): Group for which to construct the character table, with conjugacy classes in proper order
%   field ('R' or 'C'): Real or complex field
%   data (struct): JSON data of the ``realCharacterTable`` or ``complexCharacterTable`` JSON field
    characters = cellfun(@(m) m.', data.characters, 'uniform', 0); % characters is a 1D cell of 1D cell
    characters = [characters{:}].'; % which we merge into a 2D cell
    characters = replab.cyclotomic(characters); % we parse the contents
    args = cell(1, 0);
    if isfield(data, 'irreps')
        n = size(characters, 1);
        irreps = cell(1, n);
        id = data.irreps;
        for i = 1:n
            images = id{i};
            for j = 1:length(images)
                img = images{j};
                img = cellfun(@(m) m.', img, 'uniform', 0);
                img = [img{:}].';
                images{j} = replab.cyclotomic(img);
            end
            if length(images) == 0
                irreps{i} = group.trivialRep(field, 1);
            else
                irreps{i} = group.repByImages(field, size(images{1}, 1), 'images', images);
            end
        end
        args = horzcat(args, {'irreps', irreps});
    end
    if isfield(data, 'irrepNames')
        args = horzcat(args, {'irrepNames', data.irrepNames});
    end
    if isfield(data, 'classNames')
        args = horzcat(args, {'classNames', data.classNames});
    end
    if strcmp(field, 'R')
        C = replab.RealCharacterTable(group, group.conjugacyClasses, characters, args{:});
    else
        C = replab.ComplexCharacterTable(group, group.conjugacyClasses, characters, args{:});
    end
end
