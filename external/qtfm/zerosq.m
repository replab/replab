function E = zerosq(varargin)
% ZEROSQ   N-by-N quaternion matrix of zeros. Takes the same parameters as
% the Matlab function ZEROS (q.v.).

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

E = quaternion(zeros(varargin{:}));

% $Id: zerosq.m 1004 2017-11-15 17:14:09Z sangwine $

