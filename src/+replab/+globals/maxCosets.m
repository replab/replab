function value = maxCosets(newValue)
% Maximum number of cosets to enumerate during the Todd-Coxeter algorithm
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
        storedValue = 2^22;
    end
    value = storedValue;
end
