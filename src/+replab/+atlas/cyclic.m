function A = cyclic(n)
% Constructs the cyclic group of order n
    assert(n >= 2);
    name = sprintf('Cyclic group C(%d) of order %d', n, n);
    % Permutation realization
    X = [2:n 1];
    % standard presentation
    % < x | x^n = 1 >
    ag = replab.AbstractGroup({'x'}, {['x^' num2str(n)]}, 'permutationGenerators', {X}, 'order', vpi(n));
    ct = replab.ct.CyclicCharacterTable(n);
    assert(ct.group == ag.permutationGroup);
    A = replab.AtlasEntry(name, ag, ct.imap(ag.niceMorphism.inverse));
end
