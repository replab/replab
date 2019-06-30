function [Q, R] = qr(A)
% QR decomposition.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(2, 2)

% Reference:
%
% Angelika Bunse-Gerstner, Ralph Byers, Volker Mehrmann,
% 'A Quaternion QR Algorithm',
% Numerische Mathematik, 55 (1), January 1989, 83-95,
% DOI:10.1007/BF01395873
%
% Gene H. Golub and Charles van Loan, 'Matrix Computations', 3rd edition,
% Johns Hopkins University Press, 1996, section 5.1.4. ISBN 0-8018-5414-8.
% (This describes how to implement a Householder transformation without an
% explicit Householder matrix.)
%
% The code below is only loosely based on the above references. It owes
% more to the QTFM function 'bidiagonalise.m', which applies Householder
% transformations from both sides, to zero out both the lower and upper
% triangles of the matrix to be diagonalised. Here it is only necessary to
% do this from the left, to zero out the lower triangular part. As with
% bidiagonalisation, the product of the Householder matrices gives the
% unitary matrix Q. The handling of the rectangular cases is based on the
% method used in the QTFM function 'lu.m', which draws upon Golub and van
% Loan.

[m, n] = size(A);

N = min(m, n); % The number of elements on the diagonal of A.

Q = eyeq(m);   % Q is square with the same number of rows as A.
R = A;         % R is the same size as A.

for j = 1:N
    
    % Compute and apply a left Householder transformation to the jth
    % column (part of the column for the second and subsequent columns).
    
    s0 = substruct('()', {j:m, j}); % R(j:m, j)
    [h, zeta] = householder_vector(subsref(R, s0), eye(m - j + 1, 1));

    s1 = substruct('()', {j:m, j:n}); T = subsref(R, s1); % T = R(j:m, j:n)
    
    % R(j:m, j:n)     = (1./zeta) .* (T - h * (T' * h)');
    R = subsasgn(R, s1, (1./zeta) .* (T - h * (T' * h)'));
    
    s2 = substruct('()', {j:m, ':'}); T = subsref(Q, s2); % T = Q(j:m, :)
    
    % Q(j:m, :)       = (1./zeta) .* (T - h * (T' * h)');
    Q = subsasgn(Q, s2, (1./zeta) .* (T - h * (T' * h)'));
end

Q = Q'; % The algorithm above yields a Hermitian transpose, so fix it.

end

% $Id: qr.m 1004 2017-11-15 17:14:09Z sangwine $
