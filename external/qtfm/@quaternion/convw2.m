function C = convw2(A, B)
% CONVW2 2D convolution with wrap around. Unlike the standard convolution
% implemented by the CONV2 function, the result from this function matches
% that which would be obtained by a Fourier domain product, that is it
% treats the input as if extended periodically, as a tiling of the plane.

% Computes the convolution of A with B. A may be a single quaternion array,
% or a pair of quaternion arrays of the same size, in which case a
% double-sided convolution A1 * B * A2 is computed. To pass a pair of
% arrays for A, A must be a cell array with two elements (as is done for
% the conv function (q.v.)).

% Copyright (c) 2013-2015 Stephen J. Sangwine and Todd A. Ell.
% See the file : Copyright.m for further details.

narginchk(2, 2), nargoutchk(0, 1)

if ~isa(B, 'quaternion')
    error('Second input parameter must be a quaternion array')
end

if iscell(A)
    
    if length(A) ~= 2
        error('If a cell array, the first parameter must have 2 elements.')
    end
    
    L = A{1}; R = A{2}; % Extract two matrices from the cell array.
    
    if any(size(L) ~= size(R))
        error('The elements of the cell array must have the same size.')
    end
    
    if ~ismatrix(L) || ~ismatrix(R)
        error('Both elements of the cell array must be matrices')
    end
    
    if ~isa(L, 'quaternion') || ~isa(R, 'quaternion')
        error('Elements of the cell array must be quaternions')
    end
    
else
    if ~isa(A, 'quaternion') || ~isa(B, 'quaternion')
        error('Both parameters must be quaternions if the first is not a cell array')
    end

    if any(size(A) ~= size(B))
        error('The two parameters must have the same size (for now)')
    end
    
    if ~ismatrix(A) || ~ismatrix(B)
        error('Both parameters must be matrices')
    end
    
    L = A; % Set the left array to be A, and supply a default for R in the
    R = 1; % case where A is not a cell array (R is then a scalar 1).
end

[r, c] = size(L);

C = quaternion(L); % Preallocate the result array. All will be overwritten.

RB = rot90(B, 2); % Reversed copy of B (horizontally and vertically).

for row = 1:r
    for col = 1:c
        % This is a class method, therefore normal array indexing does not
        % work. See the file implementation_notes.txt, point 8.

        %C(row, col) = sum(sum(L .* circshift(R, [row, col]) .* R));
        C = subsasgn(C, substruct('()', {row,col}), ...
                        sum(sum(L .* circshift(RB, [row, col]) .* R)));
    end
end

% $Id: convw2.m 1004 2017-11-15 17:14:09Z sangwine $
