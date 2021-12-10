function value = cosetEnumerationMethod(newValue)
% Method to use for coset enumeration
%
% RepLAB can either use 'R'elator-based or 'C'oset-table based enumeration, see Chapter 5 of
% Holt, Derek, Handbook of Computational Group Theory, Chapman & Hall/CRC, 2004, pp. 149–198
%
%
% Args:
%   newValue ({'R', 'C'} or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 'C';
    end
    value = storedValue;
end
