function [V,D,W] = realeig(A)
% Real eigendecomposition
%
% The returned expressions obey ``A*V = V*D`` ``W'*A = D*W'``
    assert(isreal(A));
    [V,D,W] = eig(A,'vector');
    ind = find(imag(D) ~= 0);
    ind = ind(1:2:end);
    if ~isempty(ind)
        T = speye(length(D));
        invT = speye(length(D));
        for j = ind(:)'
            assert(D(j) == conj(D(j+1)));
            T(j:j+1,j:j+1) = [1 1; 1i -1i];
            invT(j:j+1,j:j+1) = [1 -1i; 1 1i]/2;
        end
        V = V*invT;
        D = T*diag(sparse(D))*invT;
        W = W*invT;
    else
        D = diag(sparse(D));
    end
end
