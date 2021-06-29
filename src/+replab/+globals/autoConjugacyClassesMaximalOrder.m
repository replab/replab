function value = autoConjugacyClassesMaximalOrder(newValue)
% Maximum order for which algorithms depending on conjugacy classes are automatically used
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
        storedValue = 10000;
    end
    value = storedValue;
end
