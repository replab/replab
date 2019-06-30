function write(filename, format, a)
% WRITE    Output a quaternion array to a text file.
%
% write(filename, format, a)
%
% This function writes a text file which is designed to be easily
% read into other software, such as Maple, Mathematica, etc.  The
% first line of the file contains two integer values giving the
% number of rows and columns in the file.  The quaternion values
% then follow in raster order, that is with the column index varying
% most rapidly. Each quaternion value occupies one line in the file,
% and consists of four floating-point values, separated by spaces.
%
% The format parameter may be omitted, in which case a default is
% assumed which outputs sufficient digits to represent double values
% with no loss of accuracy. Otherwise the format is a string (see the
% standard Matlab function fprintf for details).
%
% If the parameter a is a pure quaternion, only three components are
% written per quaternion.

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(2, 3), nargoutchk(0, 0) 

if nargin == 2
    
    % Only two parameters have been supplied.  The missing one must
    % be the format, since the other two cannot be omitted. Shuffle
    % the parameters, and supply the default value for the format.

    a      = format;
    format = '%+24.16e';    
end

[r, c] = size(a);

fid = fopen(filename, 'w');

fprintf(fid, '%8u %8u\n', r, c);
fprintf(fid, format, a);
fprintf(fid, '\n'); % For some reason the last \n gets omitted in the previous line.

if fclose(fid)
    error('The file was not written successfully.');
end

% $Id: write.m 1004 2017-11-15 17:14:09Z sangwine $

