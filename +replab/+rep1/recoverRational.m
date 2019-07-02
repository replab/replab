function V = recoverRational(sub)
% Replaces the basis sub.U by a basis with integer coefficients if possible
    assert(isequal(sub.field, 'R'));
    assert(isequal(replab.rep1.realType(sub), 'R'));
    P = sub.U*sub.U';
    [~, jb] = rref(P);
    U = P(:, jb);
    [num den] = rat(U);
    G = replab.rep1.lcm(unique(den(:)));
    V = [];
    limit = 1000;
    if G < limit
        V = num .* (G./den);
        for i = 2:size(V, 2)
            vi = V(:,i);
            vi1 = vi;
            for j = 1:i-1
                vj = V(:,j);
                vi1 = vi1*dot(vj,vj) - dot(vj,vi1)*vj;
                G = replab.rep1.gcd(vi1);
                vi1 = vi1/G;
                if max(abs(vi1)) > limit
                    V = [];
                    return
                end
            end
            V(:,i) = vi1;
        end
    end
end
