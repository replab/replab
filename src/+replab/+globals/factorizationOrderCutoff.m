function value = factorizationOrderCutoff(newValue)
% Maximum order under which use orbit enumeration to factorize group elements in generator words
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
        storedValue = 65536;
    end
    value = storedValue;
end
