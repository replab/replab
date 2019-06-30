function Y = expm(X)
% Matrix exponential.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

Y = overloadm(mfilename, X);

% TODO Implement a more accurate dedicated algorithm for this function. A
% possible candidate is perhaps given in one of the two articles referenced
% by the Matlab documentation page on expm (Moler and van Loan, 1978/2003;
% or Higham, 2005).

% $Id: expm.m 1004 2017-11-15 17:14:09Z sangwine $
