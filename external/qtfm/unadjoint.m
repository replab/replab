function C = unadjoint(A, F, B)
% UNADJOINT   Given an adjoint matrix, constructs the original
%             quaternion matrix.
%
% unadjoint(A) or
% unadjoint(A, 'complex')    assumes A is a complex    adjoint matrix.
% unadjoint(A, 'real')       assumes A is a real       adjoint matrix.
% unadjoint(A, 'quaternion') assumes A is a quaternion adjoint matrix.
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
% option 'quaternion'
%
% This function reverses the result of the ADJOINT function, so
% that UNADJOINT(ADJOINT(A)) == A.

% Copyright (c) 2005, 2009, 2010 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% References: see the code for ADJOINT.

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

% We have to be sure that A is an adjoint matrix, and not just some
% arbitrary real, complex or quaternion matrix. To do this we reconstruct
% the original matrix, assuming A is an adjoint, then build the adjoint
% from the reconstruction and compare it with A. If these are not equal to
% within a tolerance, A was not an adjoint. We are also able to eliminate
% some matrices as non-adjoints simply because of their dimensions.

[r, c] = size(A);

if strcmp(F, 'complex')
    
    if isa(A, 'quaternion')
        error('A quaternion matrix cannot be a complex adjoint')
    end
    
    if ~isnumeric(A)
        error(['An adjoint matrix must be numeric, given ', class(A)])
    end
    
    if block
        A = deinterleave(A, 2);
    end

    r2 = r / 2;
    c2 = c / 2;
    
    if fix(r2) * 2 ~= r || fix(c2) * 2 ~= c
        error('Matrix is incorrect size to be a complex adjoint');
    end

    A1 = A(     1 : r2,      1 : c2); % ( A1  A2 )
    A2 = A(     1 : r2, c2 + 1 : c ); % ( A3  A4 )
    A3 = A(r2 + 1 : r,       1 : c2);
    A4 = A(r2 + 1 : r,  c2 + 1 : c );
    
    C = quaternion(A1 + A4, (A1 - A4)./1i,...
                   A2 - A3, (A2 + A3)./1i)./2;

elseif strcmp(F, 'real')

    if isa(A, 'quaternion')
        error('A quaternion matrix cannot be a real adjoint')
    end
    
    if ~isnumeric(A)
        error(['An adjoint matrix must be numeric, given ', class(A)])
    end
    
    if block
        A = deinterleave(A, 4);
    end

    r4 = r / 4;
    c4 = c / 4;
    
    if fix(r4) * 4 ~= r || fix(c4) * 4 ~= c
        error('Matrix is incorrect size to be a real adjoint');
    end
   
    W = A(         1 :     r4, 1 : c4);
    X = A(    r4 + 1 : 2 * r4, 1 : c4);
    Y = A(2 * r4 + 1 : 3 * r4, 1 : c4);
    Z = A(3 * r4 + 1 :     r , 1 : c4);
    
    C = quaternion(W, X, Y, Z);
    
else % F must be 'quaternion' since we checked it above.
    
    if block
        error('''block'' is not currently supported with ''quaternion''.')
    end

    if ~isa(A, 'quaternion')
        error(['Matrix is not a quaternion matrix,',...
               ' cannot be a quaternion adjoint'])
    end
    
    if ~isreal(A)
        error(['Matrix is not a real quaternion matrix,',...
               ' cannot be a quaternion adjoint'])
    end
    
    r2 = r / 2;
    c2 = c / 2;
    
    if fix(r2) * 2 ~= r || fix(c2) * 2 ~= c
        error('Matrix is incorrect size to be a quaternion adjoint');
    end
    
    % Rebuild the original complex quaternion matrix from the top half of
    % the adjoint.
    
    C = complex(A(1:r2, 1:c2), A(1:r2, c2 + 1:c));
end

T = adjoint(C, F);
if any(size(T) ~= size(A))
    error(['Matrix is incorrect size to be a ', F, ' adjoint']);
end

% Now we verify that A was indeed an adjoint matrix to within a tolerance.
% We cannot test for exact equality, because A may be a product or the
% result of some algorithm that theoretically yields an adjoint matrix,
% but in practice yields an almost-adjoint because of rounding errors.

RT = real(T); RA = real(A); 
IT = imag(T); IA = imag(A); % In the case of real quaternions, these will be 0.

if any(any(abs(RT - RA) .* eps > abs(RA) | abs(IT - IA) .* eps > abs(IA)))
    warning('QTFM:inaccuracy', ...
            'Matrix is not (accurately) a valid adjoint matrix.');
end
end

% TODO: the above test is not ideal. It is difficult to devise a simple
% test that will work for complexified quaternions (because the modulus
% can vanish). We can define a non-vanishing quasi-modulus and use that,
% but it might not work correctly for all cases. For the moment the above
% test will do. It will correctly raise an error in simple cases of bad
% code.

function B = deinterleave(A, S)
% Internal function to alter the layout if the block option has been given.
% S is the size of the blocks representing one element. This code reverses
% the effect of the matching function interleave in adjoint.m (q.v.). It
% works by the same calculation, but it carries out the row and column
% moves in the reverse direction, applying the computed base and offset
% indices to the input rather than the output. Comparison of the code with
% the interleave code will explain this better than a comment.

[R, C] = size(A);

M = R ./ S; % The number of matrix blocks stacked vertically within A.

T = zeros(size(A), class(A)); % Initialise an intermediate array to store
                              % re-ordered rows of the input matrix.

for r = 1:R % For each row in the input matrix ...
    base    = floor((r - 1) ./ S);
    offset  = mod((r - 1), S) * M + 1;
    T(base + offset, :) = A(r, :); % Copy a row from the input matrix.
end

N = C ./ S; % The number of matrix blocks stacked horizontally within A.

% Now repeat the deinterleaving on the columns to complete the task.

B = zeros(size(A), class(A)); % Initialise the output array.

for c = 1:C % For each column in the intermediate matrix ...
    base    = floor((c - 1) ./ S);
    offset  = mod((c - 1), S) * N + 1;
    B(:, base + offset) = T(:, c); % Copy a column from the intermediate.
end
end

% $Id: unadjoint.m 1004 2017-11-15 17:14:09Z sangwine $
