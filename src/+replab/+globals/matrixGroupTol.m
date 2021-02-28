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
    persistent value
    if nargin == 1
        value = newValue;
    elseif isempty(value)
        value = 1e-10;
    end
end
