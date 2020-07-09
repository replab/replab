function value = atlas(newValue)
% Get/set the group atlas to use for group recognition
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = replab.atlas.Standard;
    end
    value = storedValue;
end
