function disp(q)
% DISP Display array.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005, 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 0)

% The argument is not checked, since this function is called by Matlab only
% if the argument is a quaternion.  There are three cases to be handled:
% empty, a pure quaternion, a full quaternion.  In the latter two cases,
% the fields may be arrays.

d = size(q);
S = blanks(5);
if isempty(q)
    if sum(d) == 0
        S = [S '[] quaternion'];
    else
        S = [S 'Empty quaternion'];
        l = length(d);
        if l == 2
            S = [S ' matrix: '];
        else
            S = [S ' array: '];
        end
        for k = 1:l
            S = [S, num2str(d(k))];
            if k == l
                break % If we have just added the last dimension, no need
                      % for another multiplication symbol.
            end
            S = [S, '-by-'];
        end
    end
elseif isscalar(q)
    S = [S char(q)];
else
    % Non-scalar case. Build up the string piece by piece, then output it
    % when complete.

    l = length(d);
    for k = 1:l
        S = [S, num2str(d(k))];
        if k == l
            break % If we have just added the last dimension, no need for
                  % another multiplication symbol.
        end
        S = [S, 'x'];
    end
    if isempty(q.w), S = [S, ' pure'];    end
    if ~isreal(q),   S = [S, ' complex']; end
                     S = [S, ' quaternion array'];
                           
    % Add some information about the components of the quaternion, unless
    % the component type is double (the default - implied).
    
    C = class(q.x);
    if ~strcmp(C, 'double')
        S = [S, ' with ', C, ' components'];
    end
end
disp(S)

% $Id: disp.m 1004 2017-11-15 17:14:09Z sangwine $
