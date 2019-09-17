function lines = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object, that fits within the given width/height limit
%
%        obj: Object to pretty print
%
%    maxRows: maximum number of rows
% maxColumns: maximum number of columns
%
% Returns a nRows x 1 cell array of strings 'lines'; the output may be cut arbitrarily at columns
% after 'maxColumns', and for rows after 'maxRows'. It is up to calling code to shape such
% overflowing output before finally printing it.  
%
%    maxRows: maximum row size; optional parameter with default value given in replab.Parameters
% maxColumns: maximum column size; optional parameter with default value given in replab.Parameters
    if nargin < 3
        maxColumns = replab.Parameters.strMaxColumns;
    end
    if nargin < 2
        maxRows = replab.Parameters.strMaxRows;
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
