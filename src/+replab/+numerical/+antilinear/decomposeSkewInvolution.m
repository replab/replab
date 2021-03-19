function [A, J, B] = decomposeSkewInvolution(F)
% Returns a decomposition of the matrix representing an antilinear involution
%
% The matrices ``A``, ``B`` obey ``B = inv(conj(A))``, with ``J = [0 eye(n); -eye(n) 0]``, and
% and ``F = A * J * B``.
%
% Args:
%   F (double(d,d)): Matrix such that ``F*conj(F) = -eye(d)``
%
% Returns
% -------
%   A: (double(d,d))
%     First matrix in the product
%   B: (double(d,d))
%     Second matrix in the product
    d = size(F, 1);
    assert(mod(d, 2) == 0, 'Dimension must be even');
    n = d/2;
    [V, D, W] = replab.numerical.antilinear.eigSkewInvolution(F);
    D = diag(D);
    D = reshape(D(1:n), 1, n);
    % U2*F*U1 = J
    % F = inv(U2) * J * inv(U1)
    % we have conj(inv(U))*J*U = eye(d)
    A = V * diag([ones(1,n) conj(1./D)]);
    B = diag([ones(1,n) D]) * conj(W);
    J = [zeros(n,n) eye(n); -eye(n) zeros(n,n)];
end
