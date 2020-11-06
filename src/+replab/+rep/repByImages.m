function rep = repByImages(group, field, dimension, preimages, images, imagesErrorBound)
% Constructs a representation from images of group generators and their inverses
%
% Args:
%   group (`+replab.FiniteGroup`): Finite group represented, must not be trivial
%   field ({'R', 'C'}): Whether the representation if real (R) or complex (C)
%   dimension (integer): Representation dimension
%   preimages (cell(1,\*) of ``group`` elements): Preimages generating the group
%   images (cell(1,\*) of double(\*,\*), may be sparse or cyclotomic): Images of the preimages
%   imagesErrorBound (double, optional): Error bound on the images
    if nargin < 6
        imagesErrorBound = [];
    end
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
        recG = cell(1, n);
        recg = cell(1, n);
        isGenPerm = true;
        o = 1;
        for i = 1:n
            [G g] = replab.perm.GeneralizedSymmetricGroup.fromMatrix(images{i});
            if isempty(G)
                isGenPerm = false;
                break
            else
                recG{i} = G;
                recg{i} = g;
                o = lcm(o, G.m);
            end
        end
        if isGenPerm
            G = replab.perm.GeneralizedSymmetricGroup(recG{1}.n, o);
            genPerms = arrayfun(@(i) recG{i}.naturalMorphism(G).imageElement(recg{i}), 1:n, 'uniform', 0);
            rep = replab.rep.RepByImages_monomial(group, field, dimension, preimages, images, G, genPerms);
        else
            rep = replab.rep.RepByImages_exact(group, field, dimension, preimages, images);
        end
    else % ~exactImages
        if isempty(imagesErrorBound)
            warning('Error bound not provided for inexact images; using adhoc error estimation');
            imagesErrorBound = 0;
            for i = 1:n
                img = images{i};
                if isa(img, 'intval')
                    img = mid(img);
                end
                eo = group.elementOrder(preimages{i});
                imagesErrorBound = max(errorBound, norm(img^eo - eye(dimension), 'fro')/eo);
            end
        end
        rep = replab.rep.RepByImages_inexact(group, field, dimension, preimages, images, imagesErrorBound);
    end
end
