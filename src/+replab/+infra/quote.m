function q = quote(str)
% Prepares a string so that it can be written in source code
%
% It wraps the string in single quotes, and replaces any inside single quotes by two consecutive single quotes.
%
% Args:
%   str (charstring): String to quote
%
% Returns:
%   charstring: Quoted string
    q = ['''' strrep(str, '''', '''''') ''''];
end
