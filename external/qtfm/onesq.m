function E = onesq(varargin)
% ONESQ   Quaternion matrix of ones. Takes the same parameters as the
% Matlab function ONES (q.v.). NB: The vector part is zero, not ones.

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

E = quaternion(ones(varargin{:}));

% $Id: onesq.m 1004 2017-11-15 17:14:09Z sangwine $

