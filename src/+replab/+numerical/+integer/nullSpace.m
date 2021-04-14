function B = nullSpace(A)
% Returns an integer basis of the null space of an integer matrix
%
% Args:
%   A (integer(m,n)): Integer matrix
%
% Returns:
%   integer(n,r): Matrix whose columns form a basis of the null space of ``A``
    [H, U] = replab.numerical.integer.hermiteNormalForm(A');
    H = H';
    U = U';
    % now H = A*U
    mask = all(H == 0, 1);
    B = U(:,mask);
end
