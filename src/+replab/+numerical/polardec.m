function [U, H] = polardec(A)
% Returns the polar decomposition of a square matrix A = U H, where U is unitary and H is semidefinite positive
%
% Args:
%   A (double(n,n)): Matrix to decompose
%
% Returns
% -------
%   U: double(n,n)
%     Unitary matrix
%   P: double(n,n)
%     Semidefinite positive matrix
    n = size(A, 1);
    assert(size(A, 2) == n, 'Matrix must be square');
    % see https://en.wikipedia.org/wiki/Polar_decomposition#Relation_to_the_SVD
    [W,S,V] = svd(A);
    U = W*V';
    % one step of Newton iteration to have more orthogonality
    N = U'*U;
    P = U*N/2;
    U = 2*U + P*N - 3*P;
    % compute SDP factor
    H = U'*A;
end
