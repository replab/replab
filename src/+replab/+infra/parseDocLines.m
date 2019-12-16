function [pos docLines] = parseDocLines(ct, pos)
% Parses zero, one or more lines of comments and returns them with the leading % stripped
    docLines = {};
    while 1
        [res line] = ct.expect(pos, '%');
        if isempty(res)
            break
        else
            docLines{1, end+1} = line(2:end);
            pos = res;
        end
    end
end
