function ct = A5CharacterTable
    G = replab.PermutationGroup.alternating(5);
    mat = {'1','1','1','1','1'
         '3','-1','0','-E(5)-E(5)^4','-E(5)^2-E(5)^3'
         '3','-1','0','-E(5)^2-E(5)^3','-E(5)-E(5)^4'
         '4','0','1','-1','-1'
         '5','1','-1','0','0'};
    chars = replab.cyclotomic.fromStrings(mat);
    conjugacyClasses = {[1 2 3 4 5] [2 1 4 3 5] [3 1 2 4 5] [5 1 2 3 4] [4 1 2 5 3]};
    classNames = {'1+1+1+1+1' '2+2+1' '3+1+1' '5a' '5b'};
    preimages = {[5 1 2 3 4] [5 3 1 2 4]};
    images1 = cellfun(@(s) replab.cyclotomic.fromStrings(s), {{'1'} {'1'}}, 'uniform', 0);
    mat = { {'-1/2*E(5)^2-1/2*E(5)^3' '-1/2*E(5)+1/2*E(5)^2+1/2*E(5)^3-1/2*E(5)^4' '1/2*E(5)+1/2*E(5)^4'
             '-1/2' '-1/2*E(5)-1/2*E(5)^4' '1/2*E(5)^2+1/2*E(5)^3'
             '-1/2*E(5)-1/2*E(5)^4' '-1/2*E(5)^2-1/2*E(5)^3' '-1/2*E(5)+1/2*E(5)^2+1/2*E(5)^3-1/2*E(5)^4'} ...
            {'-1/2*E(5)^2-1/2*E(5)^3' '-1/2*E(5)+1/2*E(5)^2+1/2*E(5)^3-1/2*E(5)^4' '1/2*E(5)+1/2*E(5)^4'
             '1' '0' '0'
             '1/2' '1/2*E(5)+1/2*E(5)^4' '-1/2*E(5)^2-1/2*E(5)^3'} };
    images2 = cellfun(@(s) replab.cyclotomic.fromStrings(s), mat, 'uniform', 0);
    mat = { {'E(5)+2*E(5)^2+2*E(5)^3+E(5)^4' 'E(5)^2+E(5)^3' '-2*E(5)^2-2*E(5)^3'
             '1' '-E(5)^2-E(5)^3' 'E(5)^2+E(5)^3'
             'E(5)+2*E(5)^2+2*E(5)^3+E(5)^4' '-1' '-E(5)-2*E(5)^2-2*E(5)^3-E(5)^4'} ...
            {'E(5)+2*E(5)^2+2*E(5)^3+E(5)^4' 'E(5)^2+E(5)^3' '-2*E(5)^2-2*E(5)^3'
             'E(5)^2+E(5)^3' 'E(5)^2+E(5)^3' '1'
             '2*E(5)^2+2*E(5)^3' 'E(5)+2*E(5)^2+2*E(5)^3+E(5)^4' '-2*E(5)-3*E(5)^2-3*E(5)^3-2*E(5)^4'} };
    images3 = cellfun(@(s) replab.cyclotomic.fromStrings(s), mat, 'uniform', 0);
    mat = { {'0' '0' 'E(3)^2' '0'
             '1' '-1' '-E(3)' 'E(3)^2'
             '-E(3)' 'E(3)' '0' 'E(3)^2'
             '0' '1' '0' '0'} ...
            {'0' '0' '1' '0'
             '-2*E(3)-E(3)^2' '2*E(3)+E(3)^2' '-2*E(3)-E(3)^2' 'E(3)^2'
             '1' '-1' '-E(3)' 'E(3)^2'
             '2' '-1' '-E(3)' '0'} };
    images4 = cellfun(@(s) replab.cyclotomic.fromStrings(s), mat, 'uniform', 0);
    mat = { [0 0 0 1 0; -1 0 0 -1 -1; 0 -1 -1 0 -1; 0 0 1 0 0; 1 1 0 0 1] ...
            [0 0 0 1 0; 1 0 0 0 0; 0 0 1 1 1; 1 1 0 0 1; -1 -1 -1 -1 -1] };
    images5 = cellfun(@(s) replab.cyclotomic.fromDoubles(s), mat, 'uniform', 0);
    irreps = {G.repByImages('C', 1, 'preimages', preimages, 'images', images1) ...
              G.repByImages('C', 3, 'preimages', preimages, 'images', images2) ...
              G.repByImages('C', 3, 'preimages', preimages, 'images', images3) ...
              G.repByImages('C', 4, 'preimages', preimages, 'images', images4) ...
              G.repByImages('C', 5, 'preimages', preimages, 'images', images5)};
    classes = replab.ConjugacyClasses(G, cellfun(@(g) G.conjugacyClass(g), conjugacyClasses, 'uniform', 0));
    ct = replab.CharacterTable(G, classes, chars, 'classNames', classNames, 'irreps', irreps);
end
