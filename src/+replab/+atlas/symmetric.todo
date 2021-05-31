function A = symmetric(n)
% Constructs the atlas entry corresponding to the symmetric group of degree n
    assert(n > 2);
    name = sprintf('Symmetric group S(%d) of degree %d', n, n);
    % Permutation realization
    S = [2:n 1];
    T = [2 1 3:n];
    % this is the presentation from page 2100 of
    % https://www.ams.org/journals/tran/2003-355-05/S0002-9947-03-03040-X/S0002-9947-03-03040-X.pdf
    relators = {['s^' num2str(n)], 't^2', ['(s*t)^' num2str(n-1)]};
    for j = 2:floor(n/2)
        relators{1,end+1} = sprintf('(t^-1 s^-%d t s^%d)^2', j, j);
    end
    ag = replab.AbstractGroup({'s' 't'}, relators, 'permutationGenerators', {S, T}, 'order', replab.util.factorial(n));
    ct = replab.sym.SymmetricGroupCharacterTable(n);
    A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
end
