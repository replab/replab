function set_octave_path
% Set the GNU Octave path to include qtfm.
%
% Copyright (c) 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% This function adds the necessary QTFM directories to the GNU Octave path.
% It is designed to be used once on installation, but to check for previous
% runs and stop without doing any damage.

if isempty(ver('octave'))
    error('Not running under GNU Octave, doing nothing.')
end

% Check that QTFM is not already on the path. This isn't a thorough check,
% as only the root directory is checked for.

if strfind(path, qtfm_root) ~= 0
    error(['QTFM root directory is already on the GNU Octave path,', ...
           ' doing nothing.'])
end

% Add QTFM directories to the path.

addpath([qtfm_root filesep '@quaternion' filesep 'private'], '-begin');
addpath([qtfm_root filesep '@quaternion'],                   '-begin');

%THERE IS A PROBLEM HERE
%As soon as @quaternion is added to the path, Octave finds functions like
%size.m there instead of in its own directory structure. Need to know how
%Octave handles special folders with @ to stop this happening, since we need
%to put the QTFM folders near the top of the path.

addpath([qtfm_root filesep 'test'],                          '-begin');
addpath( qtfm_root,                                          '-begin');

% Output diagnostics to the user.

disp('QTFM directories have been added to the path, but the path has not')
disp('been saved. To save the path use ''savepath''.')

% $Id: set_octave_path.m 1004 2017-11-15 17:14:09Z sangwine $
