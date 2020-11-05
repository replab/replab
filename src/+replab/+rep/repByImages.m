function rep = repByImages(group, field, dimension, preimages, images)
% Constructs a representation from images of group generators and their inverses
%
% Args:
%   group (`+replab.FiniteGroup`): Finite group represented, must not be trivial
%   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
%   dimension (integer): Representation dimension
%   preimages (cell(1,\*) of ``group`` elements): Preimages generating the group
%   images (cell(1,\*) of double(\*,\*), may be sparse or cyclotomic): Images of the preimages
    assert(isa(group, 'replab.FiniteGroup'));
    assert(isa(images, 'cell') && isrow(images));
    n = length(preimages);
    assert(n == length(images));
    nicePreimages = cellfun(@(g) group.niceMorphism.imageElement(g), preimages, 'uniform', 0);
    assert(length(images) == n);
    isCyclotomic = cellfun(@(m) isa(m, 'replab.cyclotomic'), images);
    isSym = cellfun(@(m) isa(m, 'sym'), images);
    if any(isSym)
        error('We do not support symbolic images. Please use replab.cyclotomic instead.');
    end
    isIntval = cellfun(@(m) isa(m, 'intval'), images);
    isDouble = cellfun(@(m) isa(m, 'double'), images);
    isSparse = cellfun(@(m) isa(m, 'double') && issparse(m), images);
    isInteger = cellfun(@(m) isreal(m) && all(all(floor(m) == m)), images);
    exactImages = [];
    if any(isIntval) % if any is intval, convert everything to intval
        for i = 1:n
            if isCyclotomic(i) || (isDouble(i) && isInteger(i))
                images{i} = intval(images{i})
            elseif isDouble(i) && ~isInteger(i)
                images{i} = midrad(images{i}, eps(images{i})*10);
            end
        end
        exactImages = false;
    elseif any(isCyclotomic) && all(isCyclotomic | isInteger) % if any is cyclotomic and rest are exact, convert to cyclo
        for i = 1:n
            if isDouble(i)
                images{i} = replab.cyclotomic.fromDoubles(images{i});
            end
        end
        exactImages = true;
    else % otherwise, convert everything to double
        for i = 1:n
            images{i} = double(images{i});
        end
        exactImages = all(isInteger);
    end
    if exactImages
        inverseImages = replab.rep.computeInverses(group, @(x,y) x*y, preimages, images);
        isUnitary = arrayfun(@(i) full(all(all(images{i} == inverseImages{i}'))), 1:n);
        rep = replab.rep.ExactRepByImages(group, field, dimension, preimages, images, inverseImages, all(isUnitary));
    else % ~exactImages
        ind = replab.mrp.inverseIndices(group, preimages);
        isUnitary = arrayfun(@(i) ind(i) > 0 && full(all(all(images{i} == images{ind(i)}'))), 1:n);
        rep = replab.rep.InexactRepByImages(group, field, dimension, preimages, images, all(isUnitary));
    end
end
