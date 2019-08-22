function str = escape(str)
% Escapes special characters in strings
%
% Namely, it replaces the newline ``\n``, the carrier return ``\r``, 
% the tab ``\t`` and the quote ``'`` by their escape sequences.
%
% Args:
%   str (char string): String to escape
%
% Returns:
%   The escaped string
    str = strrep(str, '\n', '\\n');
    str = strrep(str, '\r', '\\r');
    str = strrep(str, '\t', '\\t');
    str = strrep(str, '''', '''''');
end
