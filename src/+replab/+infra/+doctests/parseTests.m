function doctests = parseTests(lines, errFun)
% Finds and parses the doctests in documentation
%
% Args:
%   lines (row cell array of charstring): Lines to parse
%   errFun (function_handle): Error context function called before raising an error
%                             The calling convention is ``errFun(lineNumber)`` where ``message`` is
%                             a charstring and ``lineNumber`` is the line number in ``lines`` where the error was
%                             encountered. The function may, for example, print the context in which the error
%                             occured to the console.
%
% Returns:
%   row cell array of `.DocTest`: The parsed doctests, with line numbers corresponding
%                                 to the position in the ``lines`` cell array
%
% Raises:
%   An error if the parse is unsuccesful
    errId = 'replab:docTestParseError';
    if nargin < 2
        errFun = @(l) fprintf('Error in line %d\n', l);
    end
    n = length(lines);
    content = cell(1, n);
    indent = zeros(1, n);
    % Remove leading whitespace but remember identation level for each line
    for i = 1:n
        l = lines{i};
        if isempty(l)
            indent(i) = 0;
            content{i} = '';
        else
            tokens = regexp(l, '^(\s*)(.*)', 'tokens', 'once');
            if length(tokens) == 1
                % octave doesn't produce first token if string
                % doesn't start with some spaces
                indent(i) = 0;
                content{i} = strtrim(tokens{1});
            else
                indent(i) = length(tokens{1});
                content{i} = strtrim(tokens{2});
            end
        end
    end
    % Finds doctests which begin with the 'Example:' Sphinx directive
    doctests = {};
    i = 1;
    while i <= n
        if isequal(content{i}, 'Example:')
            j = i + 1;
            while j <= n && (isempty(content{j}) || indent(j) > indent(i))
                j = j + 1;
            end
            ps = replab.infra.doctests.ParseState.fromDocTestBlock(content(i+1:j-1));
            errFun1 = @(l) errFun(l + i);
            dt = replab.infra.doctests.DocTest.parse(ps, errFun1);
            if isempty(dt)
                errFun(i);
                error(errId, 'Error parsing the doctest block');
            else
                doctests{1, end+1} = dt.mapLineNumbers(@(ln) ln + i);
            end
            i = j;
        end
        i = i + 1;
    end
end
