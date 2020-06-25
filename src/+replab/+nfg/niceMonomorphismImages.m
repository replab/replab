function I = niceMonomorphismImages(group)
% Returns all the permutations that are images of the group elements under the nice monomorphism
%
% Args:
%   group (`+replab.NiceFiniteGroup`): Nice finite group
%
% Returns:
%   double(\*,\*): Matrix with each row being an image permutation
    ds = length(group.niceMonomorphismImage(group.identity)); % domain size
    c = group.elements.toCell;
    I = zeros(length(c), ds);
    for i = 1:length(c)
        I(i,:) = group.niceMonomorphismImage(c{i});
    end
end
