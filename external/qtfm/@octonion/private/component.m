function c = component(o, n)
% Component n of an octonion (n = 1:8)

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Caution: this code may return an empty array for n = 1, if the octonion
% is pure.

switch n
    case 1
        c = s(o.a);
    case 2
        c = x(o.a);
    case 3
        c = y(o.a);
    case 4
        c = z(o.a);
    case 5
        c = s(o.b);
    case 6
        c = x(o.b);
    case 7
        c = y(o.b);
    case 8
        c = z(o.b);
end

end

% $Id: component.m 1004 2017-11-15 17:14:09Z sangwine $
