function value = maxElements(newValue)
% Maximum number of elements to return in `+replab.FiniteSet.elements`
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
        storedValue = 2^20;
    end
    value = storedValue;
end
