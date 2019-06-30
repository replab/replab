function [L, U, P] = lu(A)
% LU decomposition.
% (Quaternion overloading of standard Matlab function.)

% Copyright (c) 2010, 2011, 2016 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(2, 3)

[m, n] = size(A);

N = min(m, n); % The number of elements on the diagonal of A.

IP = 1:m; % A vector of index values for permuting rows of L.

% Reference:
%
% Algorithm 3.2.1, section 3.2.6, modified along the lines of section
% 3.2.11 of:
% Gene H. Golub and Charles van Loan, 'Matrix Computations', 3rd ed,
% Johns Hopkins University Press, 1996.

% The code below handles all three cases, m > n, m == n and m < n.

for j = 1:N
    
    % Partial pivoting: place the largest diagonal element in j:N at
    % position j (taking largest to mean largest modulus) by swapping rows.
    % We use abs twice so that the modulus of a complex quaternion array A
    % yields a real result, so the code will work for biquaternion arrays.
    
    [~, k] = max(abs(abs(subsref(diag(A), substruct('()', {j:N})))));
    if k ~=1 % If k == 1, the largest element is already at A(j, j).
        % Swap rows j and j + k - 1 in A and corresponding indices in IP.
        l = j + k - 1;
        IP([j l]) = IP([l j]);
        r1 = substruct('()', {[j l], ':'});
        r2 = substruct('()', {[l j], ':'});
        A = subsasgn(A, r1, subsref(A, r2)); % A([j l], :) = A([l j], :)
    end
    
    if j == m, break, end % If true, j+1:m would be an empty range.
    s1 = substruct('()', {j,     j});
    s2 = substruct('()', {j+1:m, j});
    %A(j+1:m, j) = A(j+1:m, j) ./ A(j, j);
    A = subsasgn(A, s2, subsref(A, s2) ./ subsref(A, s1));
    
    if j == n, break, end % If true, j+1:n would be an empty range.
    s3 = substruct('()', {j+1:m, j+1:n});
    s4 = substruct('()', {j,     j+1:n});
    %A(j+1:m, j+1:n) = A(j+1:m, j+1:n) - A(j+1:m, j) * A(j, j+1:n);
    A = subsasgn(A, s3, subsref(A, s3) - subsref(A, s2) * subsref(A, s4));
end

% The algorithm above computes L and U in place, so extract them into the
% separate results demanded by the calling profile. The diagonal of L is
% implicit in the result produced above, so we have to supply the explicit
% values which are all ones.

% L has size [m, N],     where N = min(m, n), the number of elements on the
% U has size [N, n]      diagonal of A.

U = triu(subsref(A, substruct('()', {1:N, ':'})));        % U = triu(A(1:N,:));
L = tril(subsref(A, substruct('()', {':', 1:N})), -1) ... % L = tril(A(:,1:N));
    + eye(m, N, 'like', A.x);

if nargout == 3
    % The permutation matrix P is needed. At present we have a vector of
    % index values in IP, which we can use to construct P by permuting an
    % identity matrix.

    P = eye(m); P = P(IP, :);
else
    % Output parameter P is not needed, therefore we must modify L, by
    % permuting rows using the indices in IP. This is equivalent to
    % multiplication of L on the left by P.'. This ensures that A == L * U,
    % rather than P * A == L * U, which is the case if P is returned.

    [~, IX] = sort(IP); % TODO Can we avoid the sort?
    L = subsref(L, substruct('()', {IX, ':'})); % L = L(IX, :);
end

end

% Note on subscripted references: these do not work inside a class method
% (which this is). See the file 'implementation notes.txt', item 8.

% TODO Consider adding support for the additional calling profiles where
% the permutation information is returned as a vector. The above code has
% been changed to use a vector as an intermediate step towards this (it
% doesn't appear to save any time, but it will cut down on memory usage if
% P is not returned).

% $Id: lu.m 1004 2017-11-15 17:14:09Z sangwine $
