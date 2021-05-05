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
    v = obtainedValue;
    if ischar(v) && isrow(v)
        % Handling strings
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
    elseif isa(v, 'replab.Str')
        res = v.longStr(100, 100);
        res = res(:).';
        res = cellfun(@strtrim, res, 'uniform', 0);
        res = res(~cellfun(@isempty, res));
        assertEqual(res, expectedValue, message);
    elseif ischar(v) && isrow(v)
        assert(length(expectedValue) == 1);
        assertEqual(eval(expectedValue{1}), v, message);
    elseif isscalar(v)
        switch class(v)
          case 'vpi'
            assert(length(expectedValue) == 1);
            ev = expectedValue{1};
            v = strtrim(num2str(v));
            assertEqual(ev, v, message);
          case 'double'
            assert(length(expectedValue) == 1);
            assertEqual(eval(expectedValue{1}), v, message);
          case 'replab.cyclotomic'
            assert(length(expectedValue) == 1);
            ev = expectedValue{1};
            v = strtrim(num2str(v));
            assertEqual(ev, v, message);
          case 'logical'
            assert(length(expectedValue) == 1);
            assertEqual(logical(eval(expectedValue{1})), v, message);
          otherwise
            error('Unsupported');
        end
    else
        error('Unsupported');
    end
end
