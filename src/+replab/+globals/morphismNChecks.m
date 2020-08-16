function value = morphismNChecks(newValue)
% Default number of randomized checks to perform when constructing morphisms
%
% If equal to ``inf`` and for finite groups, it performs a deterministic check.
%
% See for example `+replab.FiniteGroup.morphismByImages`
%
% Args:
%   newValue (integer or ``inf`` or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 3;
    end
    value = storedValue;
end
