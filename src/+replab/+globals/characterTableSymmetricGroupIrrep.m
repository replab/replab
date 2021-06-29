function value = characterTableSymmetricGroupIrrep(newValue)
% Default form of irreducible representations of the symmetric group
%
% Args:
%   newValue ({'specht', 'seminormal', 'orthogonal'} or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 'specht';
    end
    value = storedValue;
end
