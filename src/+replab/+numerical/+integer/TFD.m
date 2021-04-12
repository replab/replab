function [T] = TFD(K,M,N,A)
% Helper for SNF
%
%  T is set equal to 0 if Ak is congruent to 0 modulo A(k,k); otherwise,
%  Tis set equal to p, where A(p,q) is an element of Ak which is not
%  divisible by A(k,k)
    T = 0;
    for I = K:M
        for J = K:N
            if mod(A(I,J),A(K,K)) ~= 0
                T = I;
                return
            end
        end
    end
end
