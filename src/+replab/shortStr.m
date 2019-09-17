function str = shortStr(obj, maxColumns)
% Returns a one-line description of the given object, that fits within the given width maxColumns
%
%        obj: Object to pretty print
%
% maxColumns: maximum column size; if the output does not fit, it may be returned cut at an arbitrary 
%             place, provided that place is *after* the last column that fits.
%             Optional parameter with default value given in replab.Parameters
%
    if nargin < 2
        maxColumns = replab.Parameters.strMaxColumns;
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
