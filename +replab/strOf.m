function s = strOf(obj, short)
    if nargin < 2
        short = true;
    end
    try
        s = obj.str(short);
    catch
        try
            s = num2str(obj);
        catch
            s = sprintf('%s instance', class(obj));
        end
    end
end
