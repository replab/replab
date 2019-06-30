function E = zerosvo(varargin)
% ZEROSVO   N-by-N pure octonion matrix of zeros. Takes the same
% parameters as the Matlab function ZEROS (q.v.).

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

Z = zeros(varargin{:});
E = octonion(Z, Z, Z, Z, Z, Z, Z);

% $Id: zerosvo.m 1004 2017-11-15 17:14:09Z sangwine $
