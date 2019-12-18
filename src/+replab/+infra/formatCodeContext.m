function res = formatCodeContext(lines, pos, nRows)
% Formats a code fragment
%
% Args:
%   lines (row cell vector of charstring): Source code lines
%   pos (integer): Line to highlight
%   nRows (integer): How many lines to display in total
%
% Returns:
%   row cell vector of charstring: Formatted code fragment
    start = pos - floor(nRows/2);
    start = max(start, 1);
    to = min(start + nRows - 1, length(lines));
    range = start:pos;
    numbers = arrayfun(@(i) sprintf('  %d: ', i), range, 'uniform', 0);
    w = find((start:to) == pos);
    if ~isempty(w)
        l = numbers{w};
        l(1) = '*';
        numbers{w} = l;
    end
    res = replab.str.align(horzcat(numbers.', lines(range).'), 'rl');
end
