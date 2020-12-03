function value = verbosity(newValue)
% Get/set the verbosity level
%
% Args:
%   newValue (integer, optional): New value; if empty queries the global variable. 0 is silent, 1 is minimal, 2 shows diagnotics/iterations
%
% Returns:
%   integer: The stored value
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = 0;
    end
    value = storedValue;
end
