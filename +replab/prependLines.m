function lines = prependLines(lines, spaces)
    if nargin < 2
        spaces = '  ';
    end
    lines = strsplit(lines, newline);
    if length(lines) == 0
        lines = '';
        return
    end
    lines{1} = [spaces lines{1}];
    lines = strjoin(lines, [newline spaces]);
end
