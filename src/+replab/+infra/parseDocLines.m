function [pos docLines docLineNumbers] = parseDocLines(ct, pos)
% Parses zero, one or more lines of comments and returns them with the leading % stripped
    docLines = {};
    docLineNumbers = [];
    while 1
        [res line] = ct.expect(pos, '%');
        if isempty(res)
            break
        else
            docLines{1, end+1} = line(2:end);
            docLineNumbers(1, end+1) = pos;
            pos = res;
        end
    end
end
