function [V, D, W] = eigSkewInvolution(J)
% Eigendecomposition of a matrix representing an antilinear skew-involution
%
% Let F be an antilinear involution described by a matrix J such that ``F(x) = J * conj(x)``.
% As we have ``F(F(x)) = -x`` (skew), we have ``F(J*conj(x)) = J*conj(J)*x = -x`` and thus
% ``J*conj(J) = conj(J)*J = -eye(d)``.
%
% Let ``v`` be an eigenvector of ``J`` with eigenvalue ``l``: ``J*v = l*v``.
% Then we have ``conj(J)*J*v = -v`` and ``conj(J)*l*v = -v``, thus ``conj(J)*v = -(1/l)*v``, which
% means ``J*conj(v) = -conj(1/l)*conj(v)``.
%
% Thus, for every eigenvector ``v`` with associated eigenvalue ``l``, there is an eigenvector
% ``conj(v)`` with associated eigenvalue ``-conj(1/l)``.
%
% As ``J`` is invertible, the eigenvalue is nonzero; also, the dimension of ``J`` must be even.
%
% We return an eigenvalue decomposition such that ``V(:,1:n)`` and ``V(:,n+1:2*n)`` are conjugate,
% and with ``DD = diag(D)``, ``DD(1:n) = 1./conj(DD(n+1:2*n))``, with ``n = d/2``.
%
% Args:
%   J (double(d,d)): Matrix such that J*conj(J) = eye(d)
%
% Returns
% -------
%   V: double(d,d)
%     Eigenvectors
%   D: double(d,d)
%     Diagonal matrix of eigenvalues
    isSkewSymmetric = all(all(J == -J.'));
    if isSkewSymmetric
        [V, D] = schur(J);
    else
        [V, D] = eig(J);
    end
    %  J*V = V*D
    D = diag(D);
    D = D(:).';
    [~, ind1] = sort(real(D));
    [~, ind2] = sort(-real(1./conj(D)));
    ind1 = ind1(:)';
    ind2 = ind2(:)';
    assert(all(ind1 ~= ind2));
    ind = ind1(ind1 < ind2);
    V = [V(:,ind) conj(V(:,ind))];
    D = [D(ind) -1./conj(D(ind))];
    D = diag(D);
    if isSkewSymmetric
        W = V';
    else
        W = inv(V);
    end
end
