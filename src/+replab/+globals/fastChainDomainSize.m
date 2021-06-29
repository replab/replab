function value = fastChainDomainSize(newValue)
% Maximum domainSize for fast chain construction
%
% See `+replab.PermutationGroup.partialChain`
%
% Args:
%   newValue (integer or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 1000;
    end
    value = storedValue;
end
