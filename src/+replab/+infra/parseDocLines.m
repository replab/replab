function [ps docLines] = parseDocLines(ps)
% Parses zero, one or more lines of comments and returns them with the leading % stripped
    docLines = {};
    while 1
        [res line] = ps.expect('COMMENT');
        if isempty(res)
            break
        else
            docLines{1, end+1} = line(2:end);
            ps = res;
        end
    end
end
