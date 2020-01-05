function str = headerStr(obj)
% Returns a tiny one-line description of the object type without inspecting its contents
%
% Args:
%   obj: Object to pretty print
%
% Returns:
%   charstring: String representation
    if isa(obj, 'replab.Str')
        try
            str = obj.headerStr;
            return
        catch
        end
    end
    str = replab.str.headerStr(obj);
end
