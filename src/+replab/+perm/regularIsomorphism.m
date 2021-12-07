function iso = regularIsomorphism(group)
% Computes the regular isomorphism for a given permutation group
%
% Args:
%   group (`+replab.PermutationGroup`): Permutation group
%
% Returns:
%   `+replab.FiniteIsomorphism`: Regular isomorphism
    o = group.order;
    assert(o < 1e6);
    o = double(o);
    images = cell(1, group.nGenerators);
    E = group.elementsSequence;
    for i = 1:group.nGenerators
        g = group.generator(i);
        img = zeros(1, o);
        for j = 1:o
            img(j) = double(E.find(group.compose(g, E.at(j))));
        end
        images{i} = img;
    end
    iso = group.isomorphismByImages(replab.S(o), 'preimages', group.generators, 'images', images);
end
