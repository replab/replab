function filtered = filterEvalcOutput(lines)
    filtered = {};
    for i = 1:length(lines)
        l = lines{i};
        l = regexprep(l, '<a +href=".*?>', '');
        l = regexprep(l, '</a>', '');
        l = regexprep(l, '.\x08', '');
        l = regexprep(l, char(215), 'x');
        l = strrep(l, 'ans =', '');
        l = strtrim(l);
        if length(l) > 0 && ~isequal(l, 'logical')
            filtered{1, end+1} = l;
        end
    end
end
