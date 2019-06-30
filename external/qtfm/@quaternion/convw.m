function C = convw(A, B)
% CONVW Convolution with wrap around. Unlike the standard convolution
% implemented by the CONV function, the result from this function matches
% that which would be obtained by a Fourier domain product, that is it
% treats the input as if extended periodically.

% Computes the convolution of A with B. A may be a single quaternion array,
% or a pair of quaternion arrays of the same size, in which case a
% double-sided convolution A1 * B * A2 is computed. To pass a pair of
% arrays for A, A must be a cell array with two elements (as is done for
% the conv function (q.v.).

% Copyright (c) 2013-2015 Stephen J. Sangwine and Todd A. Ell.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isa(A, 'quaternion') || ~isa(B, 'quaternion')
    error('Both input parameters must be quaternion arrays')
end

if iscell(A)
    
    if length(A) ~= 2
        error('If a cell array, the first parameter must have 2 elements.')
    end
    
    L = A{1}; R = A{2}; % Extract two vectors from the cell array.
    
    if any(size(L) ~= size(R))
        error('The elements of the cell array must have the same size.')
    end
    
    if ~isvector(L)
        error('Parameters inside the cell array must be vectors.')
    end

else
    if any(size(A) ~= size(B))
        error('The two parameters must have the same size, for now.')
    end
    
    if ~isvector(A)
        error('Parameters must be vectors')
    end
    
    L = A; % Set the left array to be A, and supply a default for R.
    R = ones(size(A));
end

C = quaternion(L); % Preallocate the result array. All will be overwritten.

RB = rot90(B, 2); % Reversed copy of B.

for index = 1:length(L)
    % This is a class method, therefore normal array indexing does not
    % work. See the file implementation_notes.txt, point 8.
    T = substruct('()', {index});
    if isrow(L)
        %C(index + 1) = sum(L .* circshift(R, [1, index]) .* R);
        C = subsasgn(C, T, sum(L .* circshift(RB, [1, index]) .* R));
    else
        % The vectors are column vectors.
        C = subsasgn(C, T, sum(L .* circshift(RB, [index, 1]) .* R));        
    end
end

% $Id: convw.m 1004 2017-11-15 17:14:09Z sangwine $
