function s = cellStr(obj, maxColumns, r, c)
% Returns a string representation of an element in a vector/matrix/cell array
%
% Delegates the pretty printing to replab.str.shortStr.
%
% The indexing is given by (r, c) for 2D structures; or just by (r) for vectors, then c is omitted
    if iscell(obj)
        if nargin < 4
            s = replab.str.shortStr(obj{r}, maxColumns);
        else
            s = replab.str.shortStr(obj{r, c}, maxColumns);
        end
    else
        if nargin < 4
            s = replab.str.shortStr(obj(r), maxColumns);
        else
            s = replab.str.shortStr(obj(r, c), maxColumns);
        end
    end
end
