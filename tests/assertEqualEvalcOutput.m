function assertEqualEvalcOutput(standardOutput, variables, values, expectedLines, sourceFilename, lineNumber)
% Compares the output  of evalc with the expected lines of output
%
% Args:
%   standardOutput (charstring): Output string of `evalc`
%   variables (cell(1,\*) of charstring): Names of the variables to check
%   values (cell(1,\*)): Values of those variables
%   expectedLines (cell(1,\*) of charstring): Lines of expected output from the doctest
%   sourceFilename (charstring): Name of the source code file where the doctest is
%   lineNumber (integer): Line number of the doctest
    link = sprintf('matlab: matlab.desktop.editor.openAndGoToLine(''%s'', %d)', sourceFilename, lineNumber);
    message = sprintf('Doctest failure on <a href="%s">line %d of %s</a>', link, lineNumber, sourceFilename);

    soLines = strsplit(standardOutput, '\n');
    expectedLines = cellfun(@strtrim, expectedLines, 'uniform', 0);
    soLines = cellfun(@strtrim, soLines, 'uniform', 0);
    expectedLines = expectedLines(~cellfun(@isempty, expectedLines));
    soLines = soLines(~cellfun(@isempty, soLines));
    l = 1; % line number
    n = length(soLines);
    assert(length(expectedLines) >= n, message);
    if n > 0
        assert(isequal(soLines, expectedLines(1:n)), message);
    end
    l = n + 1;
    for i = 1:length(variables)
        variable = variables{i};
        value = values{i};
        if ~strcmp(variable, 'ans')
            assertEqual(expectedLines{l}, [variable ' ='], message);
            l = l + 1;
        end
        if ischar(value) && isrow(value)
            assertEqual(eval(expectedLines{l}), value, message);
            l = l + 1;
        elseif isscalar(value)
            if isa(value, 'double')
                assertEqual(eval(expectedLines{l}), value, message);
                l = l + 1;
            elseif isa(value, 'logical')
                assertEqual(logical(eval(expectedLines{l})), value, message);
                l = l + 1;
            else
                error('Unsupported');
            end
        else
            error('Unsupported');
        end
    end
end
