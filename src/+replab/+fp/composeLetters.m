function z = composeLetters(x, y)
% Composes the two words described the given letters
%
% If both arguments are reduced, it guarantees that the result is reduced too.
%
% Args:
%   x (integer(1,\*)): Letters of the left hand side word
%   y (integer(1,\*)): Letters of the right hand side word
%
% Returns:
%   integer(1,\*): Letters of the product
    xe = length(x);
    ys = 1;
    yn = length(y);
    while xe > 0 && ys <= yn && x(xe) == -y(ys)
        xe = xe - 1;
        ys = ys + 1;
    end
    z = [x(1:xe) y(ys:end)]; % if empty, will give zeros(1, 0)
end
