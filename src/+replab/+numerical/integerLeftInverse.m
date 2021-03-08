function [Psi, ok] = integerLeftInverse(Phi)
% Computes, if it exists, the integer left inverse of an integer matrix
%
% This function uses integer arithmetic encoded in floating-point double values,
% so it will only work if the coefficients do not explode during the computation. In practice,
% this limits the use of this function to small matrices with small coefficients.
%
% This function implements the algorithm provided in:
% M.C. Camara, L. Rodman, I.M. Spitkovsky
% One sided invertibility of matrices over commutative rings, corona problems, and Toeplitz operators with matrix symbols
% in Linear Algebra and its Applications, Volume 459, 2014, pp. 58-82
% `<https://doi.org/10.1016/j.laa.2014.06.038>_`
%
% Example:
%   >>> A = replab.numerical.pairOfIntegerInverses(4, 1);
%   >>> A = A(:,1:3);
%   >>> B = replab.numerical.integerLeftInverse(A);
%   >>> all(all(B*A == eye(3)))
%       1
%
% Args:
%   Phi (integer(n,m)): Integer matrix with ``m <= n``
%
% Returns:
%   Psi: integer(m,n)
%     Left inverse such that ``Psi * Phi == eye(m)``
%   ok: logical
%     Whether a left inverse was succesfully computed
    n = size(Phi, 1);
    m = size(Phi, 2);
    assert(m <= n);
    K = nchoosek(1:n, m);
    N = size(K, 1);
    Delta = round(arrayfun(@(k) det(Phi(K(k,:),:)), (1:N)'));
    [d, X] = replab.numerical.xgcd(Delta);
    if d ~= 1
        Psi = [];
        ok = false;
        return
    end
    DeltaStar = X;
    Psi = zeros(m, n);
    for k = 1:N
        I = K(k,:);
        PhiIk = Phi(I,:);
        PhiIkStar = zeros(m, n);
        for q = 1:m
            for i = 1:m
                p = I(i);
                Delta_pq_PhiIk = det(PhiIk([1:i-1 i+1:m], [1:q-1 q+1:m]));
                PhiIkStar(q,p) = (-1)^(i+q)*Delta_pq_PhiIk;
            end
        end
        Psi = Psi + DeltaStar(k)*PhiIkStar; % no need for (-1)^q correction actually
    end
    ok = true;
end
