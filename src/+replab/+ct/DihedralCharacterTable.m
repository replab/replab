function ct = DihedralCharacterTable(n)
% Generates the character table for the dihedral group Dn
%
% From Fässler, A., Stiefel, E., &amp; Wong, B. D. (1992). Group theoretical methods and their applications. 
% Boston: Birkhäuser, 23-24.
% 
% Args:
%   n (integer): number of dihedral group
%
% Returns:
%   ct (`+replab.CharacterTable`)
    group = replab.DihedralGroup(n);
    ord = double(group.order);
    classes = group.conjugacyClasses;
    nclasses = length(classes);
    irrepExp = cell(1, nclasses);
    irrepExp{1} = {{'1'}, {'1'}};
    irrepExp{2} = {{'1'}, {'-1'}};
    base = ['E(', num2str(n), ')'];
    if even(n)
        n1D = 4;
        irrepExp{3} = {{'-1'}, {'1'}};
        irrepExp{4} = {{'-1'}, {'-1'}};
    else
        n1D = 2;
    end
    for j = 1:nclasses - n1D
        g1 = {[base,'^',num2str(j)], '0'; '0', [base,'^',num2str(-j)]};
        g2 = {'0', [base,'^',num2str(j)]; [base,'^',num2str(-j)], '0'};
        irrepExp{n1D+j} = {g1, g2};
    end
    chars = cell(nclasses);
    chars(1, :) = {'1'};
    chars(1:n1D, 1) = {'1'};
    if nclasses > n1D
        chars(n1D+1:nclasses, 1) = {'2'};
    end
    for i = 2:nclasses
        rep = classes{i}.representative;
        dk = group.composeN(group.generators{2}, rep(1) - 1);
        if isequal(dk, rep)
            k = rep(1) - 1;
            chars{2, i} = '1';
            if even(n)
                chars{3, i} = num2str((-1)^k);
                chars{4, i} = num2str((-1)^k);
            end
            for j = 1:nclasses - n1D
                g1 = [base, '^', num2str(mod(j*k, ord)), '+', base,'^', num2str(-1*mod(j*k, ord))];
                chars{n1D+j, i} = g1;
            end
        else
            k = n - rep(1);
            chars{2, i} = '-1';
            if even(n)
                chars{3, i} = num2str((-1)^k);
                chars{4, i} = num2str((-1)^(k+1));
            end
            for j = 1:nclasses - n1D
                chars{n1D+j, i} = '0';
            end
        end
    end
    ct = replab.CharacterTable.make(group, classes, [], chars, [], irrepExp);
end