function tf = isempty(q)
% ISEMPTY True for empty matrix.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% narginchk(1, 1), nargoutchk(0, 1)

% It is sufficient to check the x component, because if x is empty so must
% be the y and z components. We must not check the scalar component,
% because this is empty for a non-empty pure quaternion.
     
tf = isempty(q.x); 

% $Id: isempty.m 1004 2017-11-15 17:14:09Z sangwine $

