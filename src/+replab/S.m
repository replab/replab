function s = S(domainSize)
% Returns the symmetric group acting on a certain domain size
%
% Alias for `+replab.SymmetricGroup`
%
% Args:
%   domainSize (integer): Domain size, must be > 0
    s = replab.SymmetricGroup(domainSize);
end
