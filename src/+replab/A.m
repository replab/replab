function grp = A(n)
% Constructs the alternating group
%
% Args:
%   n (integer): Group degree
%
% Returns:
%   `+replab.PermutationGroup`: The alternating group of degree ``n``
    grp = replab.PermutationGroup.alternating(n);
end
