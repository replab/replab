function assertEqualValue(variableName, obtainedValue, expectedValue, sourceFilename, lineNumber)
% Compares the output  of evalc with the expected lines of output
%
% Args:
%   variableName (charstring): Variable being tested
%   obtainedValue: Variable value
%   expectedValue (cell(1,\*) of charstring): Text representation of the expected value
%   sourceFilename (charstring): Name of the source code file where the doctest is
%   lineNumber (integer): Line number of the doctest
    link = sprintf('matlab: matlab.desktop.editor.openAndGoToLine(''%s'', %d)', sourceFilename, lineNumber);
    message = sprintf('Doctest failure for variable %s on <a href="%s">line %d of %s</a>', variableName, link, lineNumber, sourceFilename);
    v = obtaine
    switch class(v)
      case 'char'
        assert(isrow(v), 'Only row char vectors are supported in doctests. Multiline strings contain the line feed character.');
        % Two cases: quoted or unquoted
        if length(expectedValue) == 1 && expectedValue{1}(1) == ''''
            % the string is quoted
            ev = eval(expectedValue{1});
            assertEqual(ev, v, message);
        else
            v = strsplit(v, char(10));
            v = cellfun(@strtrim, v, 'uniform', 0);
            v = v(~cellfun(@isempty, v));
            assertEqual(expectedValue, v, message);
        end
      case 'double'
        ev = eval(['[' strjoin(expectedValue, '\n') ']']);
        assertEqual(v, ev, message);
      case 'logical'
        ev = logical(eval(['[' strjoin(expectedValue, '\n') ']']));
        assertEqual(v, ev, message);
      case 'vpi'
        assert(isscalar(v), 'Only scalar vpi are supported');
        assert(length(expectedValue) == 1);
        ev = expectedValue{1};
        v = strtrim(num2str(v));
        assertEqual(ev, v, message);
      case 'replab.cyclotomic'
        assert(length(expectedValues) == size(v, 1));
        for r = 1:length(expectedValues)
            row = replab.cyclotomic(strsplit(expectedValues, ' '));
            assertEqual(row, v(r,:));
        end
      otherwise
        if isa(v, 'replab.Str')
            res = v.longStr(100, 100);
            res = res(:).';
            res = cellfun(@strtrim, res, 'uniform', 0);
res = res(~cellfun(@isempty, res));
assertEqual(res, expectedValue, message);
        else
            error('Unsupported');
        end
    end
end
