function s = merge2(s1, s2, mergeFun)
% Merges two hashmaps represented by two structs
%
% Args:
%   s1 (struct): First struct to merge
%   s2 (struct): Second struct to merge
%   mergeFun (function_handle): Function handle @(x1, x2) that merges elements present in both ``s1`` and ``s2``
    s = s1;
    f = fieldnames(s2);
    for i = 1:length(f)
        n = f{i};
        if isfield(s1, n)
            if nargin < 3 || isempty(mergeFun)
                error(sprintf('Key %s present in both structs', n));
            end
            s.(n) = mergeFun(s1.(n), s2.(n));
        else
            s.(n) = s2.(n);
        end
    end
end
