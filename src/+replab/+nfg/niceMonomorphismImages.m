function I = niceMonomorphismImages(group)
% Returns all the permutations that are images of the group elements under the nice monomorphism
%
% Args:
%   group (`+replab.NiceFiniteGroup`): Nice finite group
%
% Returns:
%   double(\*,\*): Matrix with each row being an image permutation
    ds = length(group.niceMonomorphismImage(group.identity)); % domain size
    generators = cell(size(group.generators));
    for i = 1:length(generators)
        generators{i} = group.niceMonomorphismImage(group.generators{i});
    end
    ord = double(group.order);
    perm_group = replab.PermutationGroup(ds, generators, ord);
    I = perm_group.niceGroup.chain.matrixFromElements;
end
