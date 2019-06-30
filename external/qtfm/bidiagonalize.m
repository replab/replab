function [U, B, V] = bidiagonalize(A)
% Bidiagonalize A,  such that U * A * V = B and U' * B * V' = A. B is the
% same size as A, has no vector part, and is upper or lower bidiagonal
% depending on its shape. U and V are unitary matrices.

% Copyright (c) 2005, 2008, 2012 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(3, 3)

% References:
%
% Sangwine, S. J. and Le Bihan, N.,
% Quaternion singular value decomposition based on bidiagonalization
% to a real or complex matrix using quaternion Householder transformations,
% Applied Mathematics and Computation, 182(1), 1 November 2006, 727-738, 
% DOI:10.1016/j.amc.2006.04.032.
%
% Sangwine, S. J. and Le Bihan, N.,
% Quaternion Singular Value Decomposition based on Bidiagonalization
% to a Real Matrix using Quaternion Householder Transformations,
% arXiv:math.NA/0603251, 10 March 2006, available at http://www.arxiv.org/
%
% Gene H. Golub and Charles van Loan, 'Matrix Computations', 3rd edition,
% Johns Hopkins University Press, 1996, section 5.1.4. ISBN 0-8018-5414-8.
% (This describes how to implement a Householder transformation without an
% explicit Householder matrix.)

[r, c] = size(A);

if prod([r, c]) == 1
    error('Cannot bidiagonalize a matrix of one element.');    
end

if c <= r
    [U, B, V] = internal_bidiagonalizer(A); % Upper bidiagonal result.
else
    % This requires a lower bidiagonal result. We handle this by operating
    % on the Hermitian transpose of A. The results require some swapping
    % and transposition to get the correct results.
    
    [V, B, U] = internal_bidiagonalizer(A'); U = U'; B = B.'; V = V';
end

B = check(B); % Verify the result and convert to exactly bidiagonal form
              % with a real or complex B.

% -------------------------------------------------------------------------

function [U, B, V] = internal_bidiagonalizer(A)

[r, c] = size(A); assert(c <= r);

% Iterate over the whole matrix, dealing with one column and one row on
% each pass through the loop. See Figure 1 in the Applied Mathematics and
% Computation paper cited above for a diagrammatic representation of this
% process.

F = str2func(class(A)); % F is a function handle.

U = F(eye(r)); % Initialise the three matrix results to appropriately sized
B = A;         % matrices. These are conformant so that U * B * V has the
V = F(eye(c)); % same size as A, and U and V have the same type as A.

for i = 1:c

    % Compute and apply a left Householder transformation to the ith
    % column (part of the column for the second and subsequent columns).
    
    [h, zeta] = householder_vector(B(i:end, i), eye(r - i + 1, 1));

    T = B(i:end, i:end); B(i:end, i:end) = (1./zeta) .* (T - h * (T' * h)');
    T = U(i:end,  :   ); U(i:end,  :   ) = (1./zeta) .* (T - h * (T' * h)');
    
    if i == c, return; end % On the last column, we are done, since there
                           % is no corresponding row. See Figure 1 in the
                           % Applied Mathematics and Computation paper
                           % cited above to see why this is.
    
    % Now compute and apply a right Householder transformation to the ith
    % row. In the case of the first row, this excludes the first element,
    % and in later rows, the first i elements.
    
    j = i + 1;
  
    [h, zeta] = householder_vector(B(i, j:end), eye(1, c - i)); g = conj(h);

    T = B(i:end, j:end); B(i:end, j:end) = (T - (T * h.') * g) ./zeta;
    T = V( :   , j:end); V( :   , j:end) = (T - (T * h.') * g) ./zeta;
end

% -------------------------------------------------------------------------

function R = check(B)

% Verify results, and convert the result to exactly bidiagonal form with no
% vector part.

[r, c] = size(B);

if r == 1 || c == 1
    % The matrix is degenerate (a row or column vector) and we have to deal
    % with it differently because the Matlab diag function in this case
    % constructs a matrix instead of extracting the diagonal (how clever to
    % use the same name for both ideas!).
    D = B(1); % The first element is the diagonal. There is no super-diagonal.
    O = B(2 : end); % The rest is the off-diagonal.
elseif c <= r
    D = [diag(B); diag(B, +1)];    % Extract the diagonal and super-diagonal.
    O = tril(B, -1) + triu(B, +2); % Extract the off-diagonal part.
else
    D = [diag(B); diag(B, -1)];    % Extract the diagonal and sub-diagonal.
    O = tril(B, -2) + triu(B, +1); % Extract the off-diagonal part.
end

% Find the largest on and off bidiagonal elements. We use abs twice to
% allow for the case where B is complex (this occurs when A is a
% complexified quaternion matrix).

T1 = max(max(abs(abs(O)))); % Find the largest off bidiagonal element.
T2 =     max(abs(abs(D)));  % Find the largest bidiagonal element.

% NB T2 and/or T1 could be exactly zero (example, if A was zero, or an
% identity matrix). Therefore we do not divide one by the other, but
% instead multiply by a tolerance.

tolerance = eps .* 1.0e4; % This is empirically determined.

if T1 > T2 * tolerance
    warning('QTFM:inaccuracy', ...
            'Result of bidiagonalization was not accurately diagonal.');
    disp('Information: largest on- and off-bidiagonal moduli were:');
    disp(T2);
    disp(T1);
end

R = B;

if ~isa(R, 'quaternion') && ~isa(R, 'octonion')
    
    % We skip the remaining steps if R is not quaternion or octonion, to
    % permit the code to work for real or complex arrays, which do not have
    % vector parts, and therefore do not need the remaining steps to null
    % out inexact vector parts.

    return
end

% The code below will be run only for quaternion or octonion parameters, A,
% because of the test just done and the premature return in the case or
% real or complex arrays.
    
% Verify that the diagonal elements have neglible vector parts.

T2 = max(abs(abs(s(D)))); % Largest scalar modulus of the bidiagonal result.
T1 = max(abs(abs(v(D)))); % Largest vector modulus of the bidiagonal result.

if T1 > T2 * tolerance
    warning('QTFM:inaccuracy', ...
        'Result of bidiagonalization was not accurately scalar.')
    disp('Information: largest on-diagonal scalar and vector moduli were:');
    disp(T2);
    disp(T1);
end

if r == 1 || c == 1
    R = s(R); % The diagonal has only one element, so we can just forget
    % the off-diagonal.
else
    R = s(R - O); % Subtract the off-diagonal part and take the scalar part
    % of the result.
end

% $Id: bidiagonalize.m 1004 2017-11-15 17:14:09Z sangwine $
