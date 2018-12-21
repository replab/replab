function U = findUnitary(rep1, rep2)
% Implements the algorithm of
% M. Mozrzymas, M. Studziński, and M. Horodecki, J. Phys. A: Math. Theor. 47, 505203 (2014).
    group = rep1.group;
    n = rep1.dimension;
    found = false;
    G = double(group.elements.size);
    for a = 1:n
        for b = 1:n
            sab = 0;
            for i = 1:G
                g = group.elements.at(i);
                r1 = rep1.image(g);
                r2 = rep2.image(g)'; % inverse
                sab = sab + r1(a,a)*r2(b,b);
            end
            if sab > 0
                found = true;
                break
            end
        end
        if found
            break
        end
    end
    U = zeros(n, n);
    for i = 1:n
        for j = 1:n
            for k = 1:G
                g = group.elements.at(k);
                r1 = conj(rep1.image(g));
                r2 = rep2.image(g);
                U(i,j) = U(i,j) + r1(a,i)*r2(b,j);
            end
        end
    end
    U = U / sqrt(sab) * sqrt(n/G);
end
