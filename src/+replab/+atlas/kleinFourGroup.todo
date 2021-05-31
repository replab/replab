function A = kleinFourGroup
% Constructs the atlas entry corresponding to the klein four-group
    name = sprintf('Klein four-group of order %d', 4);
    % Permutation realization
    X = [2,1,4,3];
    A = [3,4,1,2];
    % Presentation from the groupprops wiki
    % < x, a | a^2 = x^2 = 1, x a x^-1 = a^-1 >
    relators = {'a^2' 'x^2' 'x a x^-1 a'};
    ag = replab.AbstractGroup({'x' 'a'}, relators, 'permutationGenerators', {X, A}, 'order', vpi(4));
    classes = replab.ConjugacyClasses(ag, cellfun(@(g) ag.conjugacyClass(g), {'1' 'x' 'a' 'x a'}, 'uniform', 0));
    chars = replab.cyclotomic([1 1 1 1; 1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1]);
    classNames = {'1' 'x' 'a' 'xa'};
    irrepNames = {'trivial' 'ker x' 'ker a' 'ker xa'};
    irreps = {ag.repByImages('C', 1, 'images', {replab.cyclotomic(1) replab.cyclotomic(1)}) ...
              ag.repByImages('C', 1, 'images', {replab.cyclotomic(1) replab.cyclotomic(-1)}) ...
              ag.repByImages('C', 1, 'images', {replab.cyclotomic(-1) replab.cyclotomic(1)}) ...
              ag.repByImages('C', 1, 'images', {replab.cyclotomic(-1) replab.cyclotomic(-1)})};
    ct = replab.ComplexCharacterTable(ag, 'C', classes, chars, 'classNames', classNames, 'irrepNames', irrepNames, 'irreps', irreps);
    A = replab.AtlasEntry(name, ag, ct);
end
