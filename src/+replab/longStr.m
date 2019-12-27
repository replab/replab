function lines = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object, that fits within the given width/height limit
%
% Returns a ``nRows x 1`` cell array of strings ``lines``; the output may be cut arbitrarily at columns
% after ``maxColumns``, and for rows after ``maxRows``. It is up to
% calling code to shape such overflowing output before finally printing it.  
%
% Args:
%   obj: Object to pretty print
%   maxRows (integer or ``[]``): Maximum row size
%                                Optional parameter with default value given in `replab.settings.strMaxRows`
%   maxColumns (integer or ``[]``): Maximum column size
%                                   Optional parameter with default value given in `replab.settings.strMaxColumns`
    if nargin < 3 || isempty(maxColumns)
        maxColumns = replab.settings.strMaxColumns;
    end
    if nargin < 2 || isempty(maxRows)
        maxRows = replab.settings.strMaxRows;
    end
    if isa(obj, 'replab.Str')
        try
            lines = obj.longStr(maxRows, maxColumns);
            return
        catch
        end
    end
    lines = replab.str.longStr(obj, maxRows, maxColumns);
end
