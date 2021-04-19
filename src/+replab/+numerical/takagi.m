function [U, D] = takagi(A, algorithm, maxSweeps)
% Returns the Autonne-Takagi factorization of the given complex symmetric matrix
%
% It returns a unitary matrix ``U`` and a real nonnegative diagonal matrix ``D`` such that
% ``A = U.'*D*U``.

% We propose three methods:
%
% * The ``'jacobi'`` algorithm is based on Jacobi rotations, see `takagi_jacobi`.
%
% * The ``'svd'`` algorithm is based on the SVD decomposition and matrix square root,
%   see `<https://doi.org/10.1016/j.amc.2014.01.170>`_.
%
% * The ``'hybrid'`` algorithm first computes a unitary using the SVD decomposition, then
%   refines it using the Jacobi algorithm.
%
% Speed-wise, we have the following ordering for random complex symmetric matrices:
% ``'svd'`` is faster than ``'hybrid'`` which is faster than ``'jacobi'``
%
% Precision-wise, the ``'svd'`` method is the worst by an order of magnitude, whereas
% ``'jacobi'`` is similar to ``'hybrid'`` with a slight advantage for ``'hybrid'``.
%
% Args:
%   A (double(\*,\*)): Square complex symmetric matrix to decompose
%   algorithm ('svd', 'jacobi', 'hybrid', optional): Algorithm, default: hybrid
%   maxSweeps (integer, optional): Jacobi algorithm, maximal number of sweeps across the matrix, default: 50
%
% Returns
% -------
%   U: double(\*,\*)
%     Unitary matrix
%   D: double(\*,\*)
%     Diagonal matrix
    n = size(A, 1);
    assert(size(A, 2) == n, 'A must be square');
    if n == 1
        U = conj(sqrt(sign(A)));
        D = abs(A);
        return
    end
    if nargin < 3 || isempty(maxSweeps)
        maxSweeps = 50;
    end
    if nargin < 2 || isempty(algorithm)
        algorithm = 'hybrid';
    end
    assert(all(all(A == A.')), 'A must be symmetric'); % A is symmetric
    switch algorithm
      case 'svd'
        % convention in https://doi.org/10.1016/j.amc.2014.01.170
        [U,Lambda,V] = svd(A);
        % we have A = A.' = U*Lambda*V'
        Z = U'*conj(V);
        Uz = U*sqrtm(Z);
        % now we have A = Uz * Lambda * Uz.'
        % however, we use the opposite transpose convention
        U = Uz.';
        D = Lambda;
      case 'jacobi'
        [U, D] = replab.numerical.takagi_jacobi(A, [], maxSweeps);
      case 'hybrid'
        [U, ~] = replab.numerical.takagi(A, 'svd');
        % U is only approximately unitary, so we perform one step of a Newton iteration
        % see https://en.wikipedia.org/wiki/Orthogonal_matrix#Nearest_orthogonal_matrix
        N = U'*U;
        P = U*N/2;
        U = 2*U + P*(N - 3*eye(n));
        % second step, perform Jacobi sweeps
        [U, D] = replab.numerical.takagi_jacobi(A, U, maxSweeps);
    end
end
