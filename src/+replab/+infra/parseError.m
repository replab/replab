function parseError(ct, pos, message, varargin)
    disp(sprintf('In %s', replab.infra.linkOpen('%s:%d', '%s:%d', ct.filename, pos)));
    disp(' ');
    disp(replab.infra.formatCodeContext(ct.lines, pos, 7));
    errorId = 'replab:parseError';
    errorMsg = sprintf(message, varargin{:});
    disp(strjoin(lines, char(10)));
    if replab.platformIsOctave
        error(errorId, errorMsg);
    else
        throwAsCaller(MException(errorId, errorMsg));
    end
end
