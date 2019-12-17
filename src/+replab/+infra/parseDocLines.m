function [pos docLines docLineNumbers] = parseDocLines(ct, pos)
% Parses zero, one or more lines of comments and returns them with the leading % stripped
%
% Also verifies that, if present, the second comment line is empty (or contains only whitespace).
    docLines = {};
    docLineNumbers = [];
    l = 1;
    while 1
        [res line] = ct.expect(pos, '%');
        if isempty(res)
            break
        else
            content = line(2:end);
            if l == 2 && ~isempty(strtrim(content))
                replab.infra.parseError(ct, pos, 'Second documentation comment line should be empty');
            end
            docLines{1,l} = content;
            docLineNumbers(1,l) = pos;
            pos = res;
        end
        l = l + 1;
    end
end
