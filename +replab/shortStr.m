function str = shortStr(obj, maxColumns)
% Returns a one-line description of the given object fitting in the given limits 
%
% Args:
%   obj: Object to pretty print
%   maxColumns (integer as double): maximum column size; optional parameter with default value
%                                   given in Settings.m
%
% Returns:
%   A one-line string representation of the object. If the output would not fit
%   in the column limit, it may be cut at an arbitrary place, provided that place
%   is *after* the last column that fits, so that it still holds that
%   ``length(str) > maxColumns``
    if nargin < 2
        maxColumns = replab.Settings.strMaxColumns;
    end
    if isa(obj, 'replab.Str')
        try
            str = obj.shortStr(maxColumns);
            return
        catch
        end
    end
    str = replab.str.shortStr(obj, maxColumns);
end
