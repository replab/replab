function tf = endsWith(str, pattern)
% Implementation of Matlab 2016b endsWith function
%
% Does not support ``'IgnoreCase'``
    if isa(str, 'cell')
        tf = cellfun(@(x) replab.compat.endsWith(x, pattern), str);
    else
        if length(str) < length(pattern)
            tf = false;
        else
            tf = isequal(str((length(str)-length(pattern)+1):end), pattern);
        end
    end
end
