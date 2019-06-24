function s = longStr(obj)
    try
        s = obj.longStr;
    catch
        if isscalar(obj)
            s = replab.str.shortStr(obj, -1);
        elseif isrow(obj)
            s = replab.str.shortStr(obj);
        elseif iscolumn(obj)
            s = replab.str.shorty
            
    end
end
