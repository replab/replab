function A = dihedral(n)
% Constructs the atlas entry corresponding to the dihedral group of order 2*n
    assert(n > 2);
    name = sprintf('Dihedral group of order %d', 2*n);
    % Permutation realization
    X = [n:-1:1];
    A = [2:n 1];
    % Presentation from the groupprops wiki
    % < x, a | a^n = x^2 = 1, x a x^-1 = a^-1 >
    relators = {['a^' num2str(n)] 'x^2' 'x a x^-1 a'};
    ag = replab.AbstractGroup({'x', 'a'}, relators, 'permutationGenerators', {X, A}, 'order', vpi(2*n));
    ct = replab.ct.DihedralCharacterTable(n);
    assert(ct.group == ag.permutationGroup);
    A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
end
