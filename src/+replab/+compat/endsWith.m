function l = endsWith(str, pattern)
% Returns whether the given string ends with the given pattern
%
% Args:
%   str (charstring): String being tested
%   pattern (charstring): Pattern to test
%
% Returns:
%   logical: Whether ``str`` has sufficient length and it ends with ``pattern`` as a substring
    l = length(str) >= length(pattern) && isequal(str((length(str)-length(pattern)+1):end), pattern);
end
