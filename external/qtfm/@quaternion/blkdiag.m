function R = blkdiag(varargin)
% BLKDIAG Construct block diagonal matrix from input arguments.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2008, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 1)

% Check that all the arguments are quaternion arrays. We don't handle the
% case where some are not, because this can be done more easily by the
% user using the quaternion constructor to convert non-quaternion
% arguments.

C = cellfun(@class, varargin, 'UniformOutput', false);
if ~all(strcmp(C, 'quaternion'))
    error('All parameters must be quaternion arrays.')
end

% Extract the components of the vector part of each input argument. The
% results here are cell arrays of the same length as varargin.

X = cellfun(@x, varargin, 'UniformOutput', false);
Y = cellfun(@y, varargin, 'UniformOutput', false);
Z = cellfun(@z, varargin, 'UniformOutput', false);

if all(cellfun(@ispure, varargin))
    
    % None of the input arguments is a full quaternion, so we can return a
    % pure quaternion result and ignore the W or scalar component.
    
    R = quaternion(blkdiag(X{:}), blkdiag(Y{:}), blkdiag(Z{:}));
else
    
    % At least one of the input arguments is a full quaternion, so we need
    % to return a full quaternion result. We use the scalar function to
    % extract the W or scalar component, because this supplies a zero
    % result for a pure argument and thus avoids us having to check the
    % arguments individually.
    
    W = cellfun(@scalar, varargin, 'UniformOutput', false);
    
    R = quaternion(blkdiag(W{:}), blkdiag(X{:}), ...
                   blkdiag(Y{:}), blkdiag(Z{:}));
end

% $Id: blkdiag.m 1004 2017-11-15 17:14:09Z sangwine $
