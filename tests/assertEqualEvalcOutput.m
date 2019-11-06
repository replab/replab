function assertEqualEvalcOutput(obtainedString, expectedLines, message)
% expectOutput Compares output of evalc with the expected lines of output
%
% Args:
%   obtainedString (char): Output string of `evalc`
%   expectedLines (cell array of char): Expected lines as strings in a cell row vector
%   message (char) : message to be printed in case of failure
    obtainedLines = filter(strsplit(obtainedString, '\n'));
    expectedLines = filter(expectedLines);
    assertEqual(obtainedLines, expectedLines, message);
end
function filtered = filter(lines)
    filtered = {};
    for i = 1:length(lines)
        l = lines{i};
        l = regexprep(l, '<a +href=".*?>', '');
        l = regexprep(l, '</a>', '');
        l = regexprep(l, '.\x08', '');
        l = strrep(l, 'ans =', '');
        l = strtrim(l);
        if length(l) > 0 && ~isequal(l, 'logical')
            filtered{1, end+1} = l;
        end
    end
end
