function s = S(domainSize)
% Returns the symmetric group acting on a certain domain size
%
% Alias for `+replab.PermutationGroup.symmetric`
%
% Args:
%   domainSize (integer): Domain size, must be >= 0
%
% Returns:
%   `+replab.PermutationGroup`: Symmetric group
    s = replab.SymmetricGroup.make(domainSize);
end
