function [code comment] = splitComment(l)
% Splits a line around the '%' character that denotes a comment, ignoring the '%' inside single quoted strings
%
% Works also when the source line contains no such character.
%
% Args:
%   l (charstring): Source line to split
%
% Returns
% -------
%   code:
%     charstring: The code line before the '%'
%   comment:
%     charstring: The comment after the '%' (if there is no comment, this return value is empty)
    q = (l == ''''); % single quotes
    insideComment = logical(bitand(cumsum(q), 1)); % boolean mask of comments, excluding final single quotes
    commentStart = and(~insideComment, l == '%');
    pos = find(commentStart);
    if ~isempty(pos)
        code = l(1:pos(1)-1);
        comment = l(pos(1)+1:end);
    else
        code = l;
        comment = '';
    end
end
