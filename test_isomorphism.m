gens1 = {[1 2 4 3 5 6 8 7 10 9 11 12] [12 11 9 10 5 6 7 8 3 4 2 1] [5 6 7 8 11 12 10 9 4 3 1 2] [1 2 3 4 7 8 6 5 12 11 9 10]};
gens2 = {[1, 2, 3, 4, 12, 11, 9, 10, 7, 8, 6, 5] [1, 2, 4, 3, 7, 8, 5, 6, 11, 12, 9, 10] [5, 6, 7, 8, 1, 2, 3, 4, 9, 10, 11, 12] [1, 2, 4, 3, 5, 6, 8, 7, 10, 9, 11, 12] [1, 2, 4, 3, 11, 12, 9, 10, 7, 8, 5, 6]};
G1 = replab.PermutationGroup.of(gens1{:});
G2 = replab.PermutationGroup.of(gens2{:});
fm = replab.mrp.FindMorphisms(G1, G2, G2, 'epimorphisms', true);
res = fm.searchUpToConjugation;
res = res{1};
G1.isMorphismByImages(G2, 'preimages', res.preimages, 'images', res.images)
%res = G1.morphismByImages(G2, 'preimages', res{1}.preimages, 'images', res{1}.images)
preimages = {[1, 2, 4, 3, 5, 6, 8, 7, 10, 9, 11, 12], [12, 11, 9, 10, 5, 6, 7, 8, 3, 4, 2, 1], [5, 6, 7, 8, 11, 12, 10, 9, 4, 3, 1, 2], [1, 2, 3, 4, 7, 8, 6, 5, 12, 11, 9, 10]};
images = {[1, 2, 4, 3, 5, 6, 8, 7, 10, 9, 11, 12], [1, 2, 3, 4, 11, 12, 10, 9, 8, 7, 5, 6], [5, 6, 7, 8, 12, 11, 9, 10, 3, 4, 2, 1], [3, 4, 2, 1, 5, 6, 7, 8, 12, 11, 9, 10]};
ch = arrayfun(@(i) [preimages{i} images{i}+12], 1:4, 'uniform', 0)
replab.PermutationGroup.of(ch{:}).order
