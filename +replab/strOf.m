function s = strOf(obj)
    try
        s = obj.str;
    catch
        try
            s = num2str(obj);
            if size(s, 1) > 1
                t = [];
                for i = 1:size(s, 1)
                    t = [t s(i,:) char(10)];
                end
                s = t;
            end
        catch
            s = sprintf('%s instance', class(obj));
        end
    end
end
