function value = strMaxColumns(newValue)
% Get/set the maximal number of columns to display when pretty printing
%
% The default stored value is ``[]``; when the stored value is ``[]``, it returns the current terminal width;
% otherwise it returns the integer stored value.
%
% Args:
%   newValue (integer or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value, or computed terminal width
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        try
            if replab.compat.isOctave
                pair = terminal_size;
                value = pair(2);
            else
                pair = matlab.desktop.commandwindow.size;
                value = pair(1);
            end
        catch ME
            value = 80;
        end
    else
        value = storedValue;
    end
end
