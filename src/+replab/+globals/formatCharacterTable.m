function value = formatCharacterTable(newValue)
% Maximum order under which use orbit enumeration to factorize group elements in generator words
%
% Args:
%   newValue ({'gap', 'plain', or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 'plain';
    end
    value = storedValue;
end
