function value = verboseInit(newValue)
% Get/set whether the initialization is verbose
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 0;
    end
    value = storedValue;
end
