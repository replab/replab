function b = isNonZeroMatrix(A, tol)
% Returns whether the 2-norm of a matrix exceeds a given tolerance
%
% This function uses various matrix norm bounds to avoid computing the comparatively expensive 2-norm
%
% Args:
%   A (double(\*,\*)): Matrix to check the norm of
%   tol (double): Tolerance
%
% Returns:
%   logical: True if norm(A, 2) > tol
    m = size(A, 1);
    n = size(A, 2);
    if m*n == 0
        b = false;
        return
    end
    % A is a m x n matrix
    % we have the norm inequalities
    % norm(A, 'inf')/sqrt(n) <= norm(A) <= norm(A, 'inf')*sqrt(m)
    % norm(A, 1)/sqrt(m) <= norm(A) <= norm(A, 1)*sqrt(n)
    norm1 = norm(A, 1);
    normInf = norm(A, 'inf');
    % norm2 lower bound
    norm2lb = max(normInf/sqrt(n), norm1/sqrt(m));
    % norm2 upper bound
    norm2ub = min(normInf*sqrt(m), norm1*sqrt(n));
    if tol < norm2lb % we are above tolerance for sure
        b = true;
    elseif norm2ub <= tol % we are below tolerance for sure
        b = false;
    else % cannot be deduced from norm_inf heuristic, compute the 2-norm
        b = norm(A, 2) > tol;
    end
end
