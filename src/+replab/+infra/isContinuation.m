function b = isContinuation(line)
% Returns whether the given line has a continuation (``...``)
%
%
% Example:
%   >>> replab.infra.isContinuation('function res = test(a, b, ... % function')
%       ans =
%       logical
%       1
%   >>> replab.infra.isContinuation('x = 2; % sets x to 2')
%       ans =
%       logical
%       0
%
% Args:
%   line (charstring): Line to check, without newline
%
% Returns:
%   logical: Whether the line has a continuation
    [code, comment] = replab.infra.splitComment(line);
    code = strtrim(code);
    b = length(code) >= 3 && isequal(code(end-2:end), '...');
end
