function str = escape(str)
% Escapes special characters in strings, like '\n', '''', etc...
    str = strrep(str, '\n', '\\n');
    str = strrep(str, '\r', '\\r');
    str = strrep(str, '\t', '\\t');
    str = strrep(str, '''', '''''');
end
