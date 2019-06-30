function A = read(filename, format)
% READ    Input a quaternion array from a text file. The format
% of the file is defined by the matching function WRITE (q.v.).
%
% A = read(filename, format) returns an array A of either full
% or pure quaternions dependent on what was in the file.  The
% format parameter may be omitted, in which case a default is
% supplied compatible with the function write. Otherwise, the
% format parameter is as for the Matlab function scanf (q.v.).

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1) 

if nargin == 1
    
    % Only one parameter has been supplied. The missing one must
    % be the format, so supply the default value, to match the
    % default value used in the function write. Note that the
    % strings differ slightly because of the differences between
    % printf and scanf.

    format = '%24e';    
end

fid = fopen(filename, 'r');

r = fscanf(fid, ' %8u', 1);
c = fscanf(fid, ' %8u', 1);
n = r .* c;

% Now read in the data into a real matrix, which we then re-arrange
% to produce a quaternion matrix A. Note that we can deduce whether
% the file contains pure or full quaternions from the number of
% floating-point values read: it will be either 3 or 4 times the
% number of quaternions indicated by n.

[R, count] = fscanf(fid, [' ', format], inf);
fclose(fid);

inc = count / n; % This value is either 3 or 4, depending on whether
                 % the file contained pure or full quaternions.

if inc ~= 3 && inc ~= 4
    error('The number of values read from the file is incorrect.');
end

if inc == 4
    A = quaternion(transpose(reshape(R(1 : inc : count - 3), c, r)), ...
                   transpose(reshape(R(2 : inc : count - 2), c, r)), ...
                   transpose(reshape(R(3 : inc : count - 1), c, r)), ...
                   transpose(reshape(R(4 : inc : count    ), c, r)));
else
    A = quaternion(transpose(reshape(R(1 : inc : count - 2), c, r)), ...
                   transpose(reshape(R(2 : inc : count - 1), c, r)), ...
                   transpose(reshape(R(3 : inc : count    ), c, r)));
end

% $Id: read.m 1004 2017-11-15 17:14:09Z sangwine $

