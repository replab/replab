function [A, B] = decomposeInvolution(J)
% Returns a decomposition of the matrix representing an antilinear involution
%
% The matrices ``A``, ``B`` obey ``B = inv(conj(A))`` and ``J = A * B``.
%
% Args:
%   J (double(d,d)): Matrix such that ``J*conj(J) = eye(d)``
%
% Returns
% -------
%   A: (double(d,d))
%     First matrix in the product
%   B: (double(d,d))
%     Second matrix in the product
    if all(all(J == J.'))
        B = replab.numerical.takagi(J);
        A = B.';
        return
    end
    d = size(J, 1);
    [V,D,n] = replab.numerical.antilinear.eigInvolution(J);
    D = diag(D);
    s = sqrt(D(1:n));
    s0 = sqrt(D(2*n+1:end));
    s = s(:).'; % to be sure of the shape
    s0 = s0(:).';
    U = conj(V)*diag([conj(s) 1./s 1./s0])*blkdiag(kron([1 1i; 1 -1i]/sqrt(2),eye(n)), eye(d-2*n));
    % we have conj(inv(U))*J*U = eye(d)
    A = conj(U);
    B = inv(U);
end
