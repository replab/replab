function CT = PermutationCharacterTable(group)
    CT = replab.CharacterTable(group,group.conjugacyClasses,replab.sym.permutationCharTableArray(group.domainSize));
end