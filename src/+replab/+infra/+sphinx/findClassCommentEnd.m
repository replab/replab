function pos = findClassCommentEnd(lines)
    pos = 2;
    while pos <= length(lines)
        l = lines{pos};
        if isempty(l) || l(1) ~= '%'
            break
        end
        pos = pos + 1;
    end
end
