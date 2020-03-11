function value = verboseDispatch(newValue)
% Get/set whether the dispatch mechanism is verbose
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = false;
    end
    value = storedValue;
end
