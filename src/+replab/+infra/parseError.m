function parseError(ct, pos, message, varargin)
    lines = cell(1, 2);
    pos = min(ct.nLines, pos); % guard in case we go beyond
    lines{1} = ['In ' replab.infra.linkOpen('%s:%d', '%s:%d', ct.filename, pos)];
    lines{2} = '';
    errorMsg = sprintf(message, varargin{:});
    from = max(1, pos - 3);
    to = min(ct.nLines, pos + 3);
    nums1 = arrayfun(@(l) sprintf('  %d: ', l), (from:pos-1)', 'uniform', 0);
    nums2 = sprintf('* %d: ', pos);
    nums3 = arrayfun(@(l) sprintf('  %d: ', l), (pos+1:to)', 'uniform', 0);
    source = replab.str.align(horzcat(vertcat(nums1, {nums2}, nums3), ct.lines(from:to)'), 'rl');
    lines = horzcat(lines, source', {''});
    errorId = 'replab:parseError';
    disp(strjoin(lines, char(10)));
    if replab.platformIsOctave
        error(errorId, errorMsg);
    else
        throwAsCaller(MException(errorId, errorMsg));
    end
end
