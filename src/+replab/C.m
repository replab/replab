function grp = C(n)
% Describes cyclic permutations over n = "domainSize" elements
%
% Args:
%   n (integer): Cyclic group order and domain size
%
% Returns:
%   `+replab.PermutationGroup`: The cyclic group of given order/domain size
    grp = replab.PermutationGroup.cyclic(n);
end
