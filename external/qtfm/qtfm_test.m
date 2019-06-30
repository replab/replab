function qtfm_test
% Run qtfm test code.
%
% Copyright (c) 2008, 2017 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

current_dir = pwd;

cd(qtfm_root)

copyright; % Display the copyright notice.

cd('test')

test

cd(current_dir);

% This file is provided so that test code can be run from the Start menu in
% versions of Matlab prior to R2012b (version 8). In later versions, it
% provides a convenient way to run the test code from the command line.

% It would be more elegant if the test code could be called directly, but
% the problem is how to specify the location in the info.xml file, since it
% will vary according to where the toolbox has been installed. Provided the
% path is set up correctly, meaning that this file can be found on the
% path, it will run the test code, because the cd command operates relative
% to the current directory, and the qtfm_root function determines where the
% QTFM files have been installed (all is that is necessary is that
% qtfm_root.m is on the current path).

% $Id: qtfm_test.m 1004 2017-11-15 17:14:09Z sangwine $
