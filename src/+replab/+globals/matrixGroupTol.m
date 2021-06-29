function value = matrixGroupTol(newValue)
% Tolerance for equality tests in continuous matrix groups
%
% See `+replab.ClassicalCompactGroup`.
%
% Args:
%   newValue (nonnegative double or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   double: Stored integer value
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    elseif isempty(storedValue)
        storedValue = 1e-10;
    end
    value = storedValue;
end
