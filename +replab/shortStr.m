function [str overLimit] = shortStr(obj, maxColumns)
% Returns a one-line description of the given object, that fits within the given width maxColumns
%
%        obj: Object to pretty print
%
% maxColumns: maximum column size before switching to crude string description
%             optional parameter with default value 80
%
% Returns the string 'str' and a Boolean 'overLimit' that states whether the output has been shortened
% because it ran over maxColumns.
    if nargin < 2
        maxColumns = 80;
    end
    if isa(obj, 'replab.Str')
        try
            [str overLimit] = obj.shortStr(maxColumns);
            return
        catch
        end
    end
    [str overLimit] = replab.str.shortStr(obj, maxColumns);
end
