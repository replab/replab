function value = verbosity(newValue)
% Get/set the verbosity level
%
% The levels are understood as follows:
%
% - 1 is INFO, the standard log level indicating what is happening
% - 2 is DEBUG, used for information that may be needed for diagnosing issues
% - 3 is TRACE, the most fine-grained level of information, expect it to be very verbose
%
% Args:
%   newValue (integer, optional): New value; if empty queries the global variable. 0 is silent, for 1/2/3 see list above
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
