function res = formatCodeContext(source, pos, nRows)
% Formats a code fragment
%
% Args:
%   source (row cell vector of charstring, `.CodeTokens` or `.SourceElement`): Source code
%   pos (integer): Line to highlight
%   nRows (integer): How many lines to display in total
%
% Returns:
%   charstring: Formatted code fragment, with lines separated by '\n'
    if iscell(source)
        lines = source;
    elseif isa(source, 'replab.infra.CodeTokens')
        lines = source.lines;
    elseif isa(source, 'replab.infra.SourceElement')
        ct = replab.infra.CodeTokens.fromFile(source.absoluteFilename);
        lines = ct.lines;
    else
        error('Invalid argument source %s', class(source));
    end
    start = pos - floor(nRows/2);
    start = max(start, 1);
    to = min(start + nRows - 1, length(lines));
    range = start:to;
    numbers = arrayfun(@(i) sprintf('  %d: ', i), range, 'uniform', 0);
    w = find((start:to) == pos);
    if ~isempty(w)
        l = numbers{w};
        l(1) = '*';
        numbers{w} = l;
    end
    res = strjoin(replab.str.align(horzcat(numbers.', lines(range).'), 'rl'), '\n');
end
