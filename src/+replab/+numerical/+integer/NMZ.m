function [I,J,T] = NMZ(K,M,N,A)
% Smith Normal Form helper
%
%  T is set equal to 0 if Ak = 0; otherwise t is set equal to abs(A(i,j)),
%  where A(i,j) is a nonzero element which is least in absolute value.
%
    for I = K:M
        for J = K:N
            T = abs(A(I,J));
            if T ~= 0
                break
            end
        end
        if T ~= 0
            break
        end
    end
    if T == 0
        return
    end
    for U = K:M
        for V = K:N
            S = abs(A(U,V));
            if (S ~= 0) && (S < T)
                T = S;
                I = U;
                J = V;
            end
        end
    end
end
