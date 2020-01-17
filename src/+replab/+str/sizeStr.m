function s = sizeStr(sz)
% Pretty prints size information about an object
%
% Args:
%   sz: a vector of integers
%
% Returns 's', a text representation of the form '2 x 2 x 2'.
    s = strjoin(arrayfun(@(i) num2str(i), sz, 'uniform', 0'), ' x ');
end
