function [V, D, n] = eigInvolution(J)
% Eigendecomposition of a matrix representing an antilinear involution
%
% Let F be an antilinear involution described by a matrix J such that ``F(x) = J * conj(x)``.
% As we have ``F(F(x)) = x``, we have ``F(J*conj(x)) = J*conj(J)*x = x`` and thus
% ``J*conj(J) = conj(J)*J = eye(d)``.
%
% Let ``v`` be an eigenvector of ``J`` with eigenvalue ``l``: ``J*v = l*v``.
% Then we have ``conj(J)*J*v = v`` and ``conj(J)*l*v = v``, thus ``conj(J)*v = (1/l)*v``, which
% means ``J*conj(v) = conj(1/l)*conj(v)``.
%
% Thus, for every eigenvector ``v`` with associated eigenvalue ``l``, there is an eigenvector
% ``conj(v)`` with associated eigenvalue ``conj(1/l)``.
%
% For each eigenvalue, there are two cases:
%
% - The eigenvector is real, ``v == conj(v)`` and thus ``l == conj(1/l)`` so that ``abs(l) == 1``.
% - The eigenvector is complex, ``v ~= conj(v)`` and thus there is an associated pair of eigenvectors.
%
% We return an eigenvalue decomposition such that ``V(:,1:n)`` and ``V(:,n+1:2*n)`` are conjugate, and
% ``V(:,2*n+1:end)`` is real, and with ``DD = diag(D)``, ``DD(1:n) = 1./conj(DD(n+1:2*n))``,
% while ``abs(DD(2*n+1:end)) = 1``
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
%   n: integer
%     Number of conjugate eigenvectors
    [V, D] = eig(J);
    %  J*V = V*D
    D = diag(D);
    [~, ind1] = sort(real(D));
    [~, ind2] = sort(real(1./conj(D)));
    realInd = ind1(ind1 == ind2);
    conjInd1 = ind1(ind1 < ind2);
    conjInd2 = ind1(ind1 > ind2);
    realInd = realInd(:)';
    conjInd1 = conjInd1(:)';
    conjInd2 = conjInd2(:)';
    ind = [conjInd1 conjInd2 realInd];
    V = V(:,ind);
    D = D(ind);
    n = length(conjInd1);
    V(:,2*n+1:end) = real(V(:,2*n+1:end));
    V(:,n+1:2*n) = conj(V(:,1:n));
    D(n+1:2*n) = 1./conj(D(1:n));
    D = diag(D);
end
