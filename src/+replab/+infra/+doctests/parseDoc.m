function doctests = parseDoc(doc)
% Finds and parses the doctests in the documentation of an object
%
% Args:
%   doc (`replab.infra.Doc`): Documentation
%
% Returns:
%   row cell array of `.DocTest`: The parsed doctests
%
% Raises:
%   A warning if some of the parses are unsuccesful.
    n = doc.nLines;
    content = cell(1, n);
    indent = zeros(1, n);
    % Remove leading whitespace but remember identation level for each line
    for i = 1:n
        l = doc.line(i);
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
            dt = replab.infra.doctests.DocTest.parseDocTest(doc, ps, i);
            if isempty(dt)
                warning(sprintf('Error while parsing Example: block at line %d'), i);
            else
                doctests{1, end+1} = dt;
            end
            i = j;
        end
        i = i + 1;
    end
end
