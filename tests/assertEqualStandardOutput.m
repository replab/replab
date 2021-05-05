function assertEqualStandardOutput(standardOutput, expectedLines, sourceFilename, lineNumber)
% Compares the standard output of evalc with the expected lines of output
%
% Args:
%   standardOutput (charstring): Output string of `evalc`
%   expectedLines (cell(1,\*) of charstring): Lines of expected output from the doctest
%   sourceFilename (charstring): Name of the source code file where the doctest is
%   lineNumber (integer): Line number of the doctest
    link = sprintf('matlab: matlab.desktop.editor.openAndGoToLine(''%s'', %d)', sourceFilename, lineNumber);
    message = sprintf('Doctest failure on <a href="%s">line %d of %s</a>', link, lineNumber, sourceFilename);

    soLines = strsplit(standardOutput, '\n');
    soLines = cellfun(@strtrim, soLines, 'uniform', 0);
    soLines = soLines(~cellfun(@isempty, soLines));
    if ~isequal(soLines, expectedLines)
        fprintf('Mismatch in standard output in doctest.\n\n');
        fprintf('Expected output:\n%s\n', strjoin(expectedLines, char(10)));
        fprintf('Test output:\n%s\n', strjoin(soLines, char(10)));
        error(message);
    end
end
