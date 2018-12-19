function lines = prependLines(lines, spaces)
    if nargin < 2
        spaces = '  ';
    end
    lines = strsplit(lines, char(10));
    if length(lines) == 0
        lines = '';
        return
    end
    lines{1} = [spaces lines{1}];
    lines = strjoin(lines, [char(10) spaces]);
end
