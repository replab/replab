function varargout = arrayfun(varargin)
% ARRAYFUN Apply a function to each element of an array
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2009 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

help quaternion/arrayfun;
unimplemented(mfilename);

% TODO Implement the arrayfun function. Possible method: use linear
% indexing and apply the function to the array, then reshape the result
% array to the same dimensions as the input array.

% TODO Noted that the Matlab arrayfun now seems to handle user-defined
% classes, however, by making this unimplemented function unavailable, a
% test reveals that the Matlab function requires 'UniformOutput' false, and
% it then returns a cell array. So it looks like a better approach would be
% to overload this function for quaternions, and do the job here. The
% Matlab function seems to use linear indexing, so the comment above is
% along the right lines. Skeleton code below:

% A = randq(3,4,5)
% B = A(:)
% C = reshape(B, size(A))
% and C is now identical to A.   May 2015

% $Id: arrayfun.m 1004 2017-11-15 17:14:09Z sangwine $
