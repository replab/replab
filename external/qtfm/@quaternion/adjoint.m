function C = adjoint(A, F, B)
% ADJOINT   Computes an adjoint matrix of the quaternion matrix A. The
% adjoint matrix has 'equivalent' properties to the original matrix (for
% example, singular values).
%
% adjoint(A) or
% adjoint(A, 'complex')    returns a complex adjoint matrix.
% adjoint(A, 'real')       returns a real    adjoint matrix.
% adjoint(A, 'quaternion') returns a quaternion adjoint (A must be a
%                          biquaternion matrix in this case).
%
% The third parameter (which may appear in the second position if the
% second is omitted), controls the layout of the adjoint, specifically
% whether the adjoint is organised in blocks by components (scalar, x, y,
% z) or with each quaternion represented as an adjoint block.
%
% adjoint(A, 'block')      returns a complex matrix with each quaternion
%                          represented by a 2-by-2 block.
% adjoint(A, 'real', 'block') returns a real adjoint matrix in which each
%                             quaternion is represented by a real adjoint 
%                             block, and similarly for other cases.
%
% There is no opposite for 'block'. 'block' is not supported with the
% option 'quaternion'.
%
% The definition of the adjoint matrix is not unique (several permutations
% of the layout are possible).

% Copyright (c) 2005, 2009, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 3), nargoutchk(0, 1)

switch nargin
    case 1, F = 'complex'; % Supply the default first parameter value.
            block = false; % 'block' was not given, so clear block flag.
    case 2, block = strcmp(F, 'block');
            if block, F = 'complex'; end % else F will be checked below.
    case 3, block = strcmp(B, 'block');
            if ~block
            error('Third parameter value must be ''block'', if present.')
            end
    otherwise
        error('Number of input parameters is inconsistent with nargin.')
end

if ~strcmp(F, 'real') && ~strcmp(F, 'complex') && ~strcmp(F, 'quaternion')
    error(['Second parameter value must be ''real'', ''complex''',...
           ' or ''quaternion''.'])
end

if strcmp(F, 'complex')

    % Reference:
    %
    % F. Z. Zhang, Quaternions and Matrices of Quaternions,
    % Linear Algebra and its Applications, 251, January 1997, 21-57.
    % DOI:10.1016/0024-3795(95)00543-9.
    % (The definition of the complex adjoint matrix is on page 29.)

    % Zhang's paper does not consider the case where the elements of the
    % quaternion are complex, but the adjoint definition is valid in this
    % case, although the four components of the quaternion are mixed among
    % the elements of the adjoint: this means that they must be unmixed by
    % the corresponding function unadjoint (using sums and differences).

    % Extract the components of A. We use scalar() and not A.w so that we
    % get zero if A is pure.
    
    W = scalar(A); X = A.x; Y = A.y; Z = A.z;

    C = [[ W + 1i .* X, Y + 1i .* Z]; ...
         [-Y + 1i .* Z, W - 1i .* X]];
     
     if block
         C = interleave(C, 2);
     end

elseif strcmp(F, 'real')
    
    % The definition of the real adjoint matrix (although not called that)
    % is given in several sources, but only for a singleton quaternion. The
    % extension to a matrix of quaternions follows easily from Zhang's
    % example in the complex case, and is easily checked to be valid. Up to
    % release 1.6 of QTFM the result returned was different, as shown in
    % the commented code to the right below (the first row and column were
    % transposed to obtain the new layout).
    
    % References:
    %
    % B. P. Ickes, 'A New Method for Performing Digital Control System
    % Attitude Computations using Quaternions', AIAA Journal, 8(1), January
    % 1970, pp13-17, American Institute of Aeronautics and Astronautics.
    % (Eqn 8 defines the layout used below.)
    %
    % Ward, J. P., 'Quaternions and Cayley numbers', Kluwer, 1997. (p91)
    %
    % Todd A. Ell, On Systems of Linear Quaternion Functions, February 2007
    % arXiv:math/0702084v1, http://www.arxiv.org/abs/math/0702084. (Eqn 3.)

    % Extract the components of A. We use scalar() and not s() so that we
    % get zero if A is pure.
    
    W = scalar(A); X = A.x; Y = A.y; Z = A.z;

    C = [[W, -X, -Y, -Z]; ... % C = [[ W,  X,  Y,  Z]; ...  Code up to and
         [X,  W, -Z,  Y]; ... %      [-X,  W, -Z,  Y]; ...  including QTFM
         [Y,  Z,  W, -X]; ... %      [-Y,  Z,  W, -X]; ...  version 1.6.
         [Z, -Y,  X,  W]];    %      [-Z, -Y,  X,  W]];
     
     if block
         C = interleave(C, 4);
     end
     
else % F must be 'quaternion' since we checked it above.
    
    if isreal(A)
        error(['The ''quaternion'' adjoint is defined for', ...
               ' biquaternion matrices only'])
    end
    
    if block
        error('''block'' is not currently supported with ''quaternion''.')
    end

    % Reference:
    %
    % Nicolas Le Bihan, Sebastian Miron and Jerome Mars,
    % 'MUSIC Algorithm for Vector-Sensors Array using Biquaternions',
    % IEEE Transactions on Signal Processing, 55(9), September 2007,
    % 4523-4533. DOI:10.1109/TSP.2007.896067.
    %
    % The quaternion adjoint is defined in equation 17 of the above paper,
    % on page 4525.
    
    RA = real(A); IA = imag(A);
    C = [RA, IA; -IA, RA];
   %C = [RA, IA; -IA, RA]; % Alternative code, note 2.

end
end

function B = interleave(A, S)
% Internal function to alter the layout if the block option has been given.
% S is the size of the blocks representing one element. See note 3 below
% for an explanation of the layout difference achieved by this code.

% TODO It is possible to do this much more neatly using subscripted
% assignment and strides, without loops. Example code is available within
% the private functions phi/tao/omega/nu elsewhere in the toolbox.

[R, C] = size(A);

M = R ./ S; % The number of matrix blocks stacked vertically within A.

T = zeros(size(A), class(A)); % Initialise an intermediate array to store
                              % re-ordered rows of the input matrix.

for r = 1:R % For each row in the intermediate matrix ...
    base    = floor((r - 1) ./ S);
    offset  = mod((r - 1), S) * M + 1;
    T(r, :) = A(base + offset, :); % Copy a row from the input matrix.
end

N = C ./ S; % The number of matrix blocks stacked horizontally within A.

% Now repeat the interleaving on the columns to complete the task.

B = zeros(size(A), class(A)); % Initialise the output array.

for c = 1:C % For each column in the result matrix ...
    base    = floor((c - 1) ./ S);
    offset  = mod((c - 1), S) * N + 1;
    B(:, c) = T(:, base + offset); % Copy a column from the intermediate.
end
end

% Note: in the case of a biquaternion matrix, the three cases 'complex',
% 'real' and 'quaternion' do not cover all the possibilities, and the
% option strings are not accurately descriptive (since the 'real' adjoint
% actually has complex elements). Two more possibilities can be obtained by
% a double call on this function, and to make clear what the possibilities
% are, they are listed below:
%
% adjoint(b, 'complex') gives a 2x2 complex block per element of b. This is
% not simple as the elements of the adjoint are a mix of the elements of b,
% but it works.
%
% adjoint(b, 'real') gives a 4x4 complex block per element of b. In this
% case the values in the 4x4 block are copies of the four complex elements
% of b.
%
% adjoint(b, 'quaternion') gives a 2x2 quaternion block per element of b.
% This adjoint can be further processed (since the result is a quaternion
% matrix) to give two more adjoints of b, as listed below.
%
% adjoint(adjoint(b, 'quaternion'), 'complex') gives a 4x4 complex block
% per element of b, but not the same complex values as above - the real and
% imaginary parts are swapped around.
%
% adjoint(adjoint(b, 'quaternion'), 'real') gives an 8x8 real block per
% element of b, where the elements of the block are elements of b, but with
% changes of sign.
%
% The last two cases might merit further study.

% Note 2. The 'quaternion' adjoint requires further study. The alternative
% line of code at 138 above is under study - there are possible
% compatibility issues between the various different adjoints, particularly
% for biquaternions. A matching change is needed in unadjoint.m of course.

% Note 3. In the case of a quaternion matrix (i.e. not a singleton
% quaternion), it is possible to construct the adjoint in two ways, quite
% apart from the permutations of layout mentioned above. These two ways are
% as above, with the quaternion matrix inserted in blocks into the adjoint
% (with negation and so on applied to the blocks); or in an interleaved
% form in which the blocks are distributed so that each quaternion is
% inserted into the adjoint as a block matrix following the pattern given
% above. For example, a quaternion matrix [q1, q2; q3, q4] could be
% assembled into an adjoint like this:
%
% [adjoint(q1), adjoint(q2); ...
%  adjoint(q3), adjoint(q4)]
%
% and this is what is returned by the 'block' option.

% $Id: adjoint.m 1004 2017-11-15 17:14:09Z sangwine $
