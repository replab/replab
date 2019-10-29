function assertEqualEvalcOutput(obtainedString, expectedLines, message)
% expectOutput Compares output of evalc with the expected lines of output
%
% Args:
%   obtainedString (char): Output string of `evalc`
%   expectedLines (cell array of char): Expected lines as strings in a cell row vector
%   message (char) : message to be printed in case of failure
    obtainedLines = strsplit(obtainedString, '\n');
    filteredLines = {};
    for i = 1:length(obtainedLines)
        l = strtrim(obtainedLines{i});
        if length(l) > 0
            filteredLines{1, end+1} = l;
        end
    end
    assertEqual(filteredLines, expectedLines, message);
end
