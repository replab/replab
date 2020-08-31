function ct = S4CharacterTable
    G = replab.PermutationGroup.symmetric(4);
    mat = {'1' '1' '1' '1' '1'; ...
           '1' '-1' '1' '1' '-1'; ...
           '2' '0' '2' '-1' '0'; ...
           '3' '-1' '-1' '0' '1'; ...
           '3' '1' '-1' '0' '-1'};
    chars = replab.cyclotomic.fromStrings(mat);
    conjugacyClasses = {[1 2 3 4] [2 1 3 4] [2 1 4 3] [3 1 2 4] [4 1 2 3]};
    classNames = {'1+1+1+1' '2+1+1' '2+2' '3+1' '4'};
    classes = replab.ConjugacyClasses(G, cellfun(@(g) G.conjugacyClass(g), conjugacyClasses, 'uniform', 0));
    ct = replab.CharacterTable(G, classes, chars, 'classNames', classNames);
end