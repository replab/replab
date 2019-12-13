function writeDocTest(fid, filename, sourceLine, elementName, testNumber, docTest)
% Writes a doctest in a test file
%
% Args:
%   fid (integer): File descriptor of the test source code file to write
%   relativeFilename (charstring): Source filename with path, relative to the root RepLAB directory
%                                  (not including a starting path separator, and starts with 'src')
%   sourceLine (integer): Source line where the documentation comment starts
%                         (not the doctest itself, but its parent block!)
%   elementName (charstring or []): If the doctest is in a class element, name of the element; or []
%   testNumber (integer): Index of the current test
%   docTest (+replab.+infra.DocTest): Parsed doctest
    if ~isempty(elementName)
        elementName(1) = upper(elementName(1));
    end
    fprintf(fid, 'function test%s%d\n', elementName, testNumber);
    fprintf(fid, '  filename = ''%s'';\n', filename);
    for i = 1:docTest.nCommands
        command = docTest.commands{i};
        command = cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], command, 'uniform', 0);
        output = docTest.outputs{i};
        output = strjoin(cellfun(@(x) ['''' strrep(x, '''', '''''') ''''], output, 'uniform', 0), ', ');
        if length(command) == 1
            fprintf(fid, '  out = evalc(%s);\n', command{1});
        else
            fprintf(fid, '  out = evalc(strjoin(%s, char(10)));\n', strjoin(command, ', '));
        end
        lineNumber = sourceLine +docTest.lineNumbers(i) - 1;
        fprintf(fid, '  assertEqualEvalcOutput(out, {%s}, filename, %d);\n', output, lineNumber);
    end
    fprintf(fid, 'end\n\n');
end
