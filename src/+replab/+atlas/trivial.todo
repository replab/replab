function A = trivial
% Constructs the atlas entry corresponding to the trivial group
    name = 'Trivial group';
    generators = cell(1, 0);
    relators = cell(1, 0);
    % Presentation is empty
    ag = replab.AbstractGroup(generators, relators, 'permutationGenerators', {[]});
    classes = ag.conjugacyClasses;
    characters = replab.cyclotomic.eye(1);
    classNames = {'id'};
    irrepNames = {'id'};
    irreps = {ag.trivialRep('C', 1)};
    ct = replab.ComplexCharacterTable(ag, 'C', classes, characters, 'irreps', irreps, 'classNames', classNames, 'irrepNames', irrepNames);
    A = replab.AtlasEntry(name, ag, ct);
end
