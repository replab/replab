function str = shortStr(obj, maxColumns)
% Returns a one-line description of the given object, that fits within the given width maxColumns
%
% Args:
%   obj: Object to pretty print
%   maxColumns (integer or ``[]``, optional): Maximum column size; if the output does not fit, 
%                                             it may be returned cut at an arbitrary place, provided 
%                                             that place is *after* the last column that fits.
%                                             Optional parameter with default value given in `+replab.+settings.strMaxColumns`
%
% Returns:
%   charstring: String representation
    if nargin < 2 || isempty(maxColumns)
        maxColumns = replab.settings.strMaxColumns;
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
