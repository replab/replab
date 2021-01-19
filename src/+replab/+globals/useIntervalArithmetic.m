function value = useIntervalArithmetic(newValue)
% Get/set whether to use interval arithmetic in computations
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = replab.init.intlab().works;
    end
    value = storedValue;
end
