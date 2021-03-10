function ct = S5CharacterTable
    G = replab.PermutationGroup.symmetric(5);
    mat = [1 1  1  1  1  1  1
           1 -1  1  1  -1  1  -1
           4 2  0  1  -1  -1  0
           4 -2  0  1  1  -1  0
           5 1  1  -1  1  0  -1
           5 -1  1  -1  -1  0  1
           6 0  -2  0  0  1  0];
    chars = replab.cyclotomic(mat);
    conjugacyClasses = {[1 2 3 4 5] [1 2 3 5 4] [1 3 2 5 4] [1 2 4 5 3] [2 1 4 5 3] [2 3 4 5 1] [1 3 4 5 2]};
    classNames = {'1+1+1+1+1' '2+1+1+1' '2+2+1' '3+1+1' '3+2' '5' '4+1'};
    irrepNames = {'trivial' 'sign' 'std' 'std*sign' '5a' '5b' 'extsqr'};
    classes = replab.ConjugacyClasses(G, cellfun(@(g) G.conjugacyClass(g), conjugacyClasses, 'uniform', 0));
    ct = replab.CharacterTable(G, classes, chars, 'classNames', classNames, 'irrepNames', irrepNames);
end
