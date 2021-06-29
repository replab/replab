function value = verifiedArithmetic(newValue)
% Get/set whether to use verified arithmetic
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = false;
    end
    value = storedValue;
end
