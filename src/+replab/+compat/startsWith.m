function l = startsWith(str, pattern)
% Returns whether the given string starts with the given pattern
%
% Args:
%   str (charstring): String being tested
%   pattern (charstring): Pattern to test
%
% Returns:
%   logical: Whether ``str`` has sufficient length and it starts with ``pattern`` as a substring
    l = length(str) >= length(pattern) && isequal(str(1:length(pattern)), pattern);
end
