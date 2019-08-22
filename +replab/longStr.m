function lines = longStr(obj, maxRows, maxColumns)
% Returns a multiline description of the given object that fits within constraints
%
% Args:
%   obj: Object to pretty print
%   maxRows (integer as double): maximum number of rows
%                                Optional parameter with default value given in Settings.m
%   maxColumns (integer as double): maximum number of columns
%                                   Optional parameter with default value given in Settings.m
%
% Returns:
%   A ``nRows x 1`` cell array of strings ``lines``; the output may be cut arbitrarily
%   at columns after ``maxColumns``, and for rows after ``maxRows``. It is up to 
%   calling code to shape such overflowing output before finally printing it.  
    if nargin < 3
        maxColumns = replab.Settings.strMaxColumns;
    end
    if nargin < 2
        maxRows = replab.Settings.strMaxRows;
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
