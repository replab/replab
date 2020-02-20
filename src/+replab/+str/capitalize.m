function s1 = capitalize(s)
% Converts the first character to an uppercase letter (if relevant) and lowercases all other characters
    if length(s) == 0
        s1 = s;
        return
    end
    h = s(1);
    t = s(2:end);
    s1 = [upper(h) lower(t)];
end
