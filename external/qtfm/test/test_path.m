function test_path
% Test code to verify that the Matlab path includes the QTFM folder.

% Copyright (c) 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This test code is not bombproof: it is simply a cursory check to trap the
% most obvious case where the user has not added the folder to the path as
% part of the installation process. TODO Extend this code to check also that
% the position of the folder in the path is appropriate (near the start).

disp('Checking Matlab search path ...')

p = path;
if isempty(findstr(p, 'qtfm'))
    error('The QTFM folder is not on the Matlab search path.');
else
    disp('Search path is OK.')
end

% $Id: test_path.m 1004 2017-11-15 17:14:09Z sangwine $
