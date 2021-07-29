function value = useReconstruction(newValue)
% Get/set whether to use the group reconstruction when computing equivariant projections
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = false;
    end
    value = storedValue;
end
