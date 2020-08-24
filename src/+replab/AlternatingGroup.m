function grp = AlternatingGroup(n)
% Constructs the alternating group
%
% Args:
%   n (integer): Group degree
%
% Returns:
%   `+replab.PermutationGroup`: The alternating group of degree ``n``
    warning('Deprecated. Use replab.PermutationGroup.alternating(n) instead of replab.AlternatingGroup(n)');
    grp = replab.PermutationGroup.alternating(n);
end
