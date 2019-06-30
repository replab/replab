function E = eyeq(varargin)
% EYEQ   N-by-N quaternion identity matrix. Takes the same parameters as
% the Matlab function EYE (q.v.).

% Copyright (c) 2008 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

E = quaternion(eye(varargin{:}));

% $Id: eyeq.m 1004 2017-11-15 17:14:09Z sangwine $
