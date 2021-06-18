function irreps = parseIrreps(group, field, generators, data)
% Reads the irreducible representations from JSON data
%
% Args:
%   group (`.FiniteGroup`): Finite group for which we parse the irreps
%   field ('R' or 'C'): Field over which the irreps are defined
%   generators (cell(1,\*) of ``group`` elements): Generators
%   data (cell(1,nIrreps) of cell(1,nGenerators) of images): Irrep images
%
% Returns:
%   cell(1,\*) of `.Rep`: Irreducible representations
    n = length(data);
    for i = 1:n
        images = data{i};
        for j = 1:length(images)
            img = images{j};
            img = cellfun(@(m) m.', img, 'uniform', 0);
            img = [img{:}].';
            images{j} = replab.cyclotomic(img);
        end
        if length(images) == 0
            irreps{i} = group.trivialRep(field, 1);
        else
            irreps{i} = group.repByImages(field, size(images{1}, 1), 'preimages', generators, 'images', images);
        end
    end
end