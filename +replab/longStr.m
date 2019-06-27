function [lines overLimit] = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object, that fits within the given width/height limit
% (this is the fallback implementation; user  'replab.longStr'
%
%        obj: Object to pretty print
%
%    maxRows: maximum number of rows (min. 3)
% maxColumns: maximum number of columns (min. around 40)
%
% Returns a cell array of strings 'lines' and a Boolean 'overLimit' that states whether the output has been shortened
% because it ran over limit.
%
%    maxRows: maximum row size; optional parameter with default value 25
% maxColumns: maximum column size; optional parameter with default value 120
    if nargin < 3
        maxColumns = 120;
    end
    if nargin < 2
        maxRows = 25;
    end
    if isa(obj, 'replab.Str')
        try
            [lines overLimit] = obj.longStr(maxRows, maxColumns);
            return
        catch
        end
    end
    [lines overLimit] = replab.str.longStr(obj, maxRows, maxColumns);
end
