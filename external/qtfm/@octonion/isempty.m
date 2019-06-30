function tf = isempty(q)
% ISEMPTY True for empty matrix.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.
     
tf = isempty(q.a);
if tf
    assert(isempty(q.b));
end

% $Id: isempty.m 1004 2017-11-15 17:14:09Z sangwine $
