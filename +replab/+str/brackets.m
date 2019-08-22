function [l r] = brackets(obj)
% Returns the pair of braces relevant for the given object
%
% Args:
%   obj: cell or numerical array
%
% Returns
% -------
%   l: char
%     Either ``'{'`` if ``obj`` is a cell array, ``'['`` otherwise
%   r: char
%     Either ``'}'`` if ``obj`` is a cell array, ``']'`` otherwise
    if iscell(obj)
        l = '{'; r = '}';
    else
        l = '['; r = ']';
    end
end
