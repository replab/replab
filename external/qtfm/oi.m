function o = oi
% oi is one of the seven octonion operators.
% oi is usually denoted by i, but this symbol is used in Matlab to represent
% the complex operator (also represented by j).

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

o = octonion(1, 0, 0, 0, 0, 0, 0);

% Implementation note:  i cannot be overloaded because it is a built-in Matlab
% operator. I cannot be used because Matlab does not distinguish between upper
% and lower case.

% $Id: oi.m 1004 2017-11-15 17:14:09Z sangwine $
