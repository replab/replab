function rep = repByImages(group, field, dimension, varargin)
% Constructs a representation from images of group generators and their inverses
%
% Arguments are the same as the `+replab.FiniteGroup.repByImages` function call.
    assert(isa(group, 'replab.FiniteGroup'), 'The group must be finite for the repByImages construction');
    assert(ismember(field, {'R' 'C'}), 'Field must be R or C for real or complex respectively.');
    assert(isscalar(dimension) && dimension == round(dimension), 'Dimension must be an integer');
    args = struct('preimages', {group.generators}, 'images', {{}}, 'imagesErrorBound', {[]});
    [args repKWargs] = replab.util.populateStruct(args, varargin);
    images = args.images;
    preimages = args.preimages;
    imagesErrorBound = args.imagesErrorBound;
    assert(isa(preimages, 'cell') && (isempty(preimages) || isrow(preimages)), 'Preimages are given as a cell row vector');
    assert(isa(images, 'cell') && (isempty(images) || isrow(images)), 'Images are given as a cell row vector');
    if isempty(images) && ~group.isTrivial
        error('Images must be provided for a non-trivial group');
    end
    assert(length(preimages) == length(images), 'Number of images does not match the number of preimages');

    % Detect inexact representations
    isDouble = cellfun(@(m) isa(m, 'double'), images);
    isInteger = zeros(1, length(images));
    for i = 1:length(images)
        img = images{i};
        switch class(img)
          case 'double'
            isInteger(i) = full(all(all(round(img) == img)));
          case 'replab.cyclotomic'
            isInteger(i) = all(all(img.isWhole));
        end
    end
    if any(isDouble & ~isInteger) || (~isempty(imagesErrorBound) && any(imagesErrorBound > 0))
        rep = replab.rep.RepByImages_inexact(group, field, dimension, preimages, images, imagesErrorBound, repKWargs{:});
        return
    end

    % Detect monomial representations
    n = length(preimages);
    GSG = cell(1, n);
    GSG_element = cell(1, n);
    isGenPerm = true;
    cycloOrder = 1;
    for i = 1:n
        I = images{i};
        % detect generalized permutation matrices
        [G g] = replab.perm.GeneralizedSymmetricGroup.fromMatrix(I);
        if isempty(G)
            isGenPerm = false;
            break
        else
            GSG{i} = G;
            GSG_element{i} = g;
            cycloOrder = lcm(cycloOrder, G.m);
        end
    end
    if isGenPerm
        G = replab.perm.GeneralizedSymmetricGroup(dimension, cycloOrder);
        genPerms = arrayfun(@(i) GSG{i}.naturalMorphism(G).imageElement(GSG_element{i}), 1:n, 'uniform', 0);
        rep = replab.rep.RepByImages_monomial(group, field, dimension, preimages, images, G, genPerms, repKWargs{:});
        return
    end
    rep = replab.rep.RepByImages_exact(group, field, dimension, preimages, images, repKWargs{:});
end
