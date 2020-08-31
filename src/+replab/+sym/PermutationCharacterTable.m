function CT = PermutationCharacterTable(group)
    cellArr = replab.sym.IntegerPartitions(group.domainSize).conjugacyClasses;
    conjClasses= replab.ConjugacyClasses(group,cellArr);
    CT = replab.CharacterTable(group,conjClasses,...
    replab.sym.permutationCharTableArray(group.domainSize));
end