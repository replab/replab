function ct = DihedralCharacterTable(n)
% Generates the character table for the dihedral group Dn
%
% From Fässler, A., Stiefel, E., & Wong, B. D. (1992). Group theoretical methods and their applications.
% Boston: Birkhäuser, 23-24.
%
% Args:
%   n (integer): number of dihedral group
%
% Returns:
%   ct (`+replab.CharacterTable`)
    group = replab.DihedralGroup(n);
    ord = 2*n;

    if even(n)
        nclasses = 4 + n/2 - 1;
    else
        nclasses = 2 + (n - 1) / 2;
    end

    % We have the group < d,s | d^n = s^2 = e, dxd^-1 = s^-1 >

    % Classes will be listed with first rotations and then reflections
    classreps = cell(1, nclasses);
    d = [2:n, 1];
    s = fliplr(1:n);
    if even(n)
        stop = n/2;
    else
        stop = (n-1)/2;
    end
    rep = 1:n;
    for i = 1:stop + 1
        classreps{i} = rep;
        rep = rep(d);
    end
    if n > 2
        classreps{i + 1} = s;
        if even(n)
            classreps{i + 2} = d(s);
        end
    end
    classarray = cellfun(@(r) group.conjugacyClass(r), classreps, 'UniformOutput', false);
    classes = replab.ConjugacyClasses(group, classarray);

    % Irreps are generated first in 1D and then in 2D
    irreps = cell(1, nclasses);
    w = replab.cyclotomic.E(n);
    irreps{1} = group.repByImages('R', 1, 'images', {1, 1});
    irreps{2} = group.repByImages('R', 1, 'images', {-1, 1});
    n1D = 2;
    if even(n)
        irreps{3} = group.repByImages('R', 1, 'images', {1, -1});
        irreps{4} = group.repByImages('R', 1, 'images', {-1, -1});
        stop = stop - 1;
        n1D = 4;
    end
    for j = 1:stop
        g1 = replab.cyclotomic.zeros(2, 2);
        g1(1, 1) = w^j;
        g1(2, 2) = w^(-j);
        g2 = replab.cyclotomic.zeros(2, 2);
        g2(1, 2) = w^j;
        g2(2, 1) = w^(-j);
        irreps{n1D + j} = group.repByImages('C', 2, 'images', {g2, g1});
    end

    % Characters can be assigned to rotations then reflections
    chars = replab.cyclotomic.zeros(nclasses, nclasses);
    chars(1, :) = replab.cyclotomic.fromDoubles(1);
    chars(1:n1D, 1) = replab.cyclotomic.fromDoubles(1);
    if nclasses > n1D
        chars(n1D+1:nclasses, 1) = replab.cyclotomic.fromDoubles(2);
    end
    if even(n)
        for k = 1:nclasses - 3
            chars(2, k+1) = replab.cyclotomic.fromDoubles(1);
            chars(3:4, k+1) = replab.cyclotomic.fromDoubles((-1)^k);
            for j = 1:stop
                pow = mod(j*k, ord);
                chars(n1D+j, k+1) = w^pow + w^(-pow);
            end
        end
        for k = nclasses - 2:nclasses - 1
            chars(2, k+1) = replab.cyclotomic.fromDoubles(-1);
            chars(3, k+1) = replab.cyclotomic.fromDoubles((-1)^k);
            chars(4, k+1) = replab.cyclotomic.fromDoubles((-1)^(k+1));
            for j = 1:stop
                chars(n1D+j, k+1) = replab.cyclotomic.fromDoubles(0);
            end
        end
    else
        for k = 1:nclasses - 2
            chars(2, k+1) = replab.cyclotomic.fromDoubles(1);
            for j = 1:stop
                pow = mod(j*k, ord);
                chars(n1D+j, k+1) = w^pow + w^(-pow);
            end
        end
        chars(2, nclasses) = replab.cyclotomic.fromDoubles(-1);
        for j = 1:stop
            chars(n1D+j, nclasses) = replab.cyclotomic.fromDoubles(0);
        end
    end

    ct = replab.CharacterTable(group, classes, irreps, chars);

end