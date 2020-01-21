function parseWarning(ct, pos, message, varargin)
    res = warning('query', 'replab:parseWarning');
    if isequal(res.state, 'on')
        disp(sprintf('In %s', replab.infra.repl.linkOpen('%s:%d', '%s:%d', ct.filename, pos)));
        disp(' ');
        disp(replab.infra.formatCodeContext(ct.lines, pos, 7));
    end
    warningId = 'replab:parseWarning';
    warning(warningId, message, varargin{:});
end
