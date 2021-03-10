function ct = A4CharacterTable
    G = replab.PermutationGroup.alternating(4);
    chars = replab.cyclotomic({'1' '1' '1' '1'; '1' '1' 'E(3)^2' 'E(3)'; '1' '1' 'E(3)' 'E(3)^2'; '3' '-1' '0' '0'});
    conjugacyClasses = {[1 2 3 4] [2 1 4 3] [3 1 2 4] [4 1 3 2]};
    classNames = {'1+1+1+1' '2+2' '3+1a' '3+1b'};
    preimages = {[4 1 3 2] [4 2 1 3]};
    images1 = cellfun(@(s) replab.cyclotomic(s), {{'1'} {'1'}}, 'uniform', 0);
    images2 = cellfun(@(s) replab.cyclotomic(s), {{'E(3)'} {'E(3)^2'}}, 'uniform', 0);
    images3 = cellfun(@(s) replab.cyclotomic(s), {{'E(3)^2'} {'E(3)'}}, 'uniform', 0);
    images4 = {[0 1 0; 0 0 -1; -1 0 0] [0 0 1; -1 0 0; 0 -1 0]};
    irreps = {G.repByImages('C', 1, 'preimages', preimages, 'images', images1) ...
              G.repByImages('C', 1, 'preimages', preimages, 'images', images2) ...
              G.repByImages('C', 1, 'preimages', preimages, 'images', images3) ...
              G.repByImages('C', 3, 'preimages', preimages, 'images', images4)};
    classes = replab.ConjugacyClasses(G, cellfun(@(g) G.conjugacyClass(g), conjugacyClasses, 'uniform', 0));
    ct = replab.CharacterTable(G, classes, chars, 'classNames', classNames, 'irreps', irreps);
end
