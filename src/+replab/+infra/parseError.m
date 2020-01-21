function parseError(ct, pos, message, varargin)
% Displays context information on the standard output and raises a parse error
%
% Args:
%   ct (`+replab.+infra.CodeTokens` or `+replab.+infra.SourceElement`): Context
%   pos (integer): Line number
%   message (charstring): Error message
%   varargin: Extra arguments for string formatting
    disp(sprintf('In %s', replab.infra.repl.linkOpen('%s:%d', '%s:%d', ct.filename, pos)));
    disp(' ');
    disp(replab.infra.formatCodeContext(ct.lines, pos, 7));
    errorId = 'replab:parseError';
    if replab.compat.isOctave
        error(errorId, message, varargin{:});
    else
        throwAsCaller(MException(errorId, message, varargin{:}));
    end
end
