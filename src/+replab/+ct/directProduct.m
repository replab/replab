function ct = directProduct(group, field)
    n = group.nFactors;
    classes = group.conjugacyClasses;
    tables = cellfun(@(f) f.characterTable(field), group.factors, 'uniform', 0);
    values = replab.cyclotomic(1);
    for i = 1:n
        factor_ct = group.factor(i).characterTable(field);
        val = factor_ct.values;
        cjc_ct = factor_ct.classes;
        grp_ct = group.factor(i).conjugacyClasses;
        values = kron(values, val(:, cjc_ct.indicesOfClasses(grp_ct)));
    end
    hasIrreps = all(cellfun(@(c) c.hasIrreps, tables));
    if ~hasIrreps
        if field == 'R'
            ct = replab.RealCharacterTable(group, classes, values);
        else
            ct = replab.ComplexCharacterTable(group, classes, values);
        end
    else
        irreps_factors = replab.util.cartesian(cellfun(@(c) c.irreps, tables, 'uniform', 0));
        nIrreps = length(irreps_factors);
        irreps = cell(1, nIrreps);
        for i = 1:nIrreps
            irrep_factors = irreps_factors{i};
            dims = cellfun(@(ir) ir.dimension, irrep_factors);
            images = cell(1, group.nGenerators);
            for j = 1:group.nGenerators
                g = group.generator(j);
                img = replab.cyclotomic(1);
                for k = 1:n
                    img = kron(img, irrep_factors{k}.image(g{k}, 'exact'));
                end
                images{j} = img;
            end
            irreps{i} = group.repByImages(field, prod(dims), 'preimages', group.generators, 'images', images, 'isIrreducible', true);
        end
        if field == 'R'
            ct = replab.RealCharacterTable(group, classes, values, 'irreps', irreps);
        else
            ct = replab.ComplexCharacterTable(group, classes, values, 'irreps', irreps);
        end
    end
end

function ct = directProduct1(ct1, ct2)
% Computes the character table of a direct product of two groups
%
% Args:
%   ct1 (`+replab.ComplexCharacterTable`): First character table
%   ct2 (`+replab.ComplexCharacterTable`): Second character table
%
% Returns:
%   `+replab.ComplexCharacterTable`: Character table of the product ``ct1.group.directProduct(ct2.group)``
    assert(isa(ct1, 'replab.ComplexCharacterTable'));
    assert(isa(ct2, 'replab.ComplexCharacterTable'));
    assert(ct1.field == 'C' && ct2.field == 'C');
    new_group = ct1.group.directProduct(ct2.group);
    % New characters are kronecker product of character matrices
    A = ct1.characters;
    B = ct2.characters;
    new_chars = replab.cyclotomic.zeros(size(A, 1) * size(B, 1), size(A, 2) * size(B, 2));
    for i = 0:size(A, 1)-1
        for j = 0:size(A, 2)-1
            new_chars(i * size(B,1) + 1:(i+1) * size(B,1), j * size(B,2) + 1:(j+1) * size(B,2)) = A(i+1, j+1) * B;
        end
    end
    % New conjugacy classes are cartesian product of input
    c1 = cellfun(@(c) c.representative, ct1.classes.classes, 'UniformOutput', false);
    c2 = cellfun(@(c) c.representative, ct2.classes.classes, 'UniformOutput', false);
    new_classes = cell(1, length(c1) * length(c2));
    for i = 0:length(c1) - 1
        new_classes(i*length(c2)+1:(i+1)*length(c2)) = cellfun(@(x) {c1{i+1}, x}, c2, 'UniformOutput', false);
    end
    classarray = cellfun(@(r) new_group.conjugacyClass(r), new_classes, 'UniformOutput', false);
    new_classes = replab.ConjugacyClasses(new_group, classarray);
    % New irreps are direct products of input irreps
    new_irreps = cell(1, length(ct1.irreps) * length(ct2.irreps));
    for i = 1:length(ct1.irreps)
        for j = 1:length(ct2.irreps)
            preimages1 = ct1.irreps{i}.preimages;
            preimages2 = ct2.irreps{j}.preimages;
            gens = cell(1, 0);
            for i = 1:length(preimages1)
                gens{1, end+1} = {preimages1{i} ct2.group.identity};
            end
            for i = 1:length(preimages2)
                gens{1, end+1} = {ct1.group.identity preimages2{i}};
            end
            B = ct2.irreps{j}.image(ct2.group.identity);
            new_images1 = cell(1, ct1.group.nGenerators);
            for k = 1:ct1.group.nGenerators
                A = ct1.irreps{i}.images_internal{k};
                % Kronecker product of cyclotomics
                new_image = replab.cyclotomic.zeros(size(A,1) * size(B, 1), size(A, 2) * size(B, 2));
                for m = 0:size(A, 1)-1
                    for n = 0:size(A, 2)-1
                        new_image(m*size(B,1)+1:(m+1)*size(B,1), n*size(B,2)+1:(n+1)*size(B,2)) = A(m+1, n+1) * B;
                    end
                end
                new_images1{k} = new_image;
            end
            new_images2 = cell(1, ct2.group.nGenerators);
            A = ct1.irreps{i}.image(ct1.group.identity);
            for k = 1:ct2.group.nGenerators
                B = ct2.irreps{j}.images_internal{k};
                % Kronecker product of cyclotomics
                new_image = replab.cyclotomic.zeros(size(A,1) * size(B, 1), size(A, 2) * size(B, 2));
                for m = 0:size(A, 1)-1
                    for n = 0:size(A, 2)-1
                        new_image(m*size(B,1)+1:(m+1)*size(B,1), n*size(B,2)+1:(n+1)*size(B,2)) = A(m+1, n+1) * B;
                    end
                end
                new_images2{k} = new_image;
            end
            new_dim = double(ct1.irreps{i}.dimension) * double(ct2.irreps{j}.dimension);
            new_irrep = new_group.repByImages('C', new_dim, 'preimages', gens, 'images', [new_images1, new_images2]);
            new_irreps{j + (i-1)*length(ct2.irreps)} = new_irrep;
        end
    end
    ct = replab.ComplexCharacterTable(new_group, 'C', new_classes, new_chars, 'irreps', new_irreps);
end
