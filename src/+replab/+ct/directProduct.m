function ct = directProduct(ct1, ct2)
% Computes the character table of a direct product of two groups
%
% Args:
%   ct1 (`+replab.CharacterTable`): First character table
%   ct2 (`+replab.CharacterTable`): Second character table
%
% Returns:
%   `+replab.CharacterTable`: Character table of the product ``ct1.group.directProduct(ct2.group)``
    assert(isa(ct1, 'replab.CharacterTable'));
    assert(isa(ct2, 'replab.CharacterTable'));
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
    ct = replab.CharacterTable(new_group, 'C', new_classes, new_chars, 'irreps', new_irreps);
end
