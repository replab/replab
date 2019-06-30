function q = q1
% Q1 Return quaternion unity (i.e. 1 + 0i + 0j + 0k).

% Copyright (c) 2005, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Up to version 2.4 of QTFM this function returned the same value as the qi
% function. From version 2.5 onwards it was altered to return unity, for
% consistency with qi, qj and qk, so that q1 + qi + qj + qk returns the
% quaternion value 1 + 1i + 1j + 1k.

q = quaternion(1, 0, 0, 0);

% $Id: q1.m 1004 2017-11-15 17:14:09Z sangwine $
