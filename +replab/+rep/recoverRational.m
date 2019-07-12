function V = recoverRational(sub)
% Replaces the basis sub.U by a basis with integer coefficients if possible
    V = [];
    P = sub.U'*sub.U;
    if replab.isNonZeroMatrix(imag(P), replab.Settings.doubleEigTol)
        return
    end
    P = real(P);
    [~, jb] = rref(P);
    U = P(:, jb);
    [num den] = rat(U);
    G = replab.rep.veclcm(unique(den(:)));
    limit = 1000;
    if G < limit
        V = num .* (G./den);
        V(:,1) = V(:,1)/replab.rep.vecgcd(V(:,1));
        for i = 2:size(V, 2)
            vi = V(:,i);
            vi1 = vi;
            for j = 1:i-1
                vj = V(:,j);
                vi1 = vi1*dot(vj,vj) - dot(vj,vi1)*vj;
                G = replab.rep.vecgcd(vi1);
                vi1 = vi1/G;
                if max(abs(vi1)) > limit
                    return
                end
            end
            V(:,i) = vi1;
        end
        V = V';
    end
end
