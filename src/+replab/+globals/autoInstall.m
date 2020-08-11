function value = autoInstall(newValue)
% Get/set whether dependencies should be automatically downloaded and installed
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = false;
    end
    value = storedValue;
end
