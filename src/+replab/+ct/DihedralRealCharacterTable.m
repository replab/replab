function ct = DihedralRealCharacterTable(n)
% Generates the character table for the dihedral group Dn
%
% From
%   https://groupprops.subwiki.org/wiki/Linear_representation_theory_of_dihedral_groups
% and
%   Fässler, A., Stiefel, E., & Wong, B. D. (1992). Group theoretical methods and their applications. Boston: Birkhäuser, 23-24.
%
% Args:
%   n (integer): number of dihedral group
%
% Returns:
%   ct (`+replab.CharacterTable`)
    group = replab.PermutationGroup.dihedral(n);
    ord = 2*n;

    if even(n)
        nclasses = 4 + n/2 - 1;
    else
        nclasses = 2 + (n - 1) / 2;
    end

    % We have the group < d,s | d^n = s^2 = e, dxd^-1 = s^-1 >

    % Classes will be listed with first rotations and then reflections
    if n == 2
        % Use the order d^0, d^1, sd^0, sd^1 to be consistent
        classreps = {[1,2,3,4], [3,4,1,2], [2,1,4,3], [4,3,2,1]};
        stop = 1;
    else
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
                classreps{i + 2} = s(d);
            end
        end
    end
    classarray = cellfun(@(r) group.conjugacyClass(r), classreps, 'UniformOutput', false);
    classes = replab.ConjugacyClasses(group, classarray);

    % Irreps are generated first in 1D and then in 2D
    irreps = cell(1, nclasses);
    irreps{1} = group.repByImages('C', 1, 'images', {1, 1});
    irreps{2} = group.repByImages('C', 1, 'images', {-1, 1});
    n1D = 2;
    if even(n)
        irreps{3} = group.repByImages('C', 1, 'images', {1, -1});
        irreps{4} = group.repByImages('C', 1, 'images', {-1, -1});
        stop = stop - 1;
        n1D = 4;
    end
    for j = 1:stop
        c = (w^j + w^(-j))/2;  % cos(2*pi*j/n)
        s = (w^j - w^(-j))/2i; % sin(2*pi*j/n)
        g1 = [c -s
              s c];
        g2 = replab.cyclotomic.zeros(2, 2);
        g2(1, 1) = 1;
        g2(2, 2) = -1;
        irreps{n1D + j} = group.repByImages('C', 2, 'images', {g2, g1});
    end

    % Characters can be assigned to rotations then reflections
    chars = replab.cyclotomic.zeros(nclasses, nclasses);
    chars(1, :) = replab.cyclotomic(1);
    chars(1:n1D, 1) = replab.cyclotomic(1);
    if nclasses > n1D
        chars(n1D+1:nclasses, 1) = replab.cyclotomic(2);
    end
    if even(n)
        for k = 1:nclasses - 3
            chars(2, k+1) = replab.cyclotomic(1);
            chars(3:4, k+1) = replab.cyclotomic((-1)^k);
            for j = 1:stop
                pow = mod(j*k, ord);
                chars(n1D+j, k+1) = w^pow + w^(-pow);
            end
        end
        for k = 0:1
            chars(2, nclasses - 1 + k) = replab.cyclotomic(-1);
            chars(3, nclasses - 1 + k) = replab.cyclotomic((-1)^k);
            chars(4, nclasses - 1 + k) = replab.cyclotomic((-1)^(k+1));
            for j = 1:stop
                chars(n1D+j, nclasses - 1 + k) = replab.cyclotomic(0);
            end
        end
    else
        for k = 1:nclasses - 2
            chars(2, k+1) = replab.cyclotomic(1);
            for j = 1:stop
                pow = mod(j*k, ord);
                chars(n1D+j, k+1) = w^pow + w^(-pow);
            end
        end
        chars(2, nclasses) = replab.cyclotomic(-1);
        for j = 1:stop
            chars(n1D+j, nclasses) = replab.cyclotomic(0);
        end
    end
    ct = replab.CharacterTable(group, classes, chars, 'irreps', irreps);
end
