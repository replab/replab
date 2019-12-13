function assertEqualEvalcOutput(obtainedString, expectedLines, sourceFilename, lineNumber)
% expectOutput Compares output of evalc with the expected lines of output
%
% Args:
%   obtainedString (char): Output string of `evalc`
%   expectedLines (cell array of char): Expected lines as strings in a cell row vector
%   message (char) : message to be printed in case of failure
    obtainedLines = filterEvalcOutput(strsplit(obtainedString, '\n'));
    expectedLines = filterEvalcOutput(expectedLines);
    link = sprintf('matlab: matlab.desktop.editor.openAndGoToLine(''%s'', %d)', sourceFilename, lineNumber);
    message = sprintf('Doctest failure on <a href="%s">line %d of %s</a>', link, lineNumber, sourceFilename);
    assertEqual(obtainedLines, expectedLines, message);
end
