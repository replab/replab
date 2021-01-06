function value = yolo(newValue)
% Get/set whether to error on unsound error bounds
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = true;
    end
    value = storedValue;
end
