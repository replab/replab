function n = numel(q, varargin)
% NUMEL   Number of elements in an array or subscripted array expression.
% (Quaternion overloading of standard Matlab function.)
%
% Copyright (c) 2005, 2008, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if isempty(q)
    n = 0; return % Return zero for an empty quaternion array.
end

if nargin == 1
    n = builtin('numel', q.x);
else
    % We have more than one argument. This means varargin represents a
    % list of index values. This function should probably never be
    % called with this parameter profile (Matlab should call the built-in
    % numel function instead), so we trap any call that does occur.

    error('Quaternion numel function called with varargin, unexpected.');

    % If we do have to handle this case, here is how it could be done:
    % n = numel(q.x, varargin);
end

% The history of this function is complicated. Prior to May 2008 it
% returned prod(size(q)), but it was changed in connection with the way
% Matlab implemented subscripted indexing, in order to support subsref and
% subsasgn, and returned 1 for any quaternion array. From Matlab R2015b, a
% new function numArgumentsFromSubscript was introduced, and by overloading
% this, we can now provide a better way for subsref and subsasgn to work,
% leaving numel to do the obvious and return the number of quaternions in
% an array.

% $Id: numel.m 1004 2017-11-15 17:14:09Z sangwine $
