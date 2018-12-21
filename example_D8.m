a = [2 3 4 1];
x = [3 2 1 4];
D8 = replab.PermutationGroup.fromGenerators(a, x);
irrFaith = D8.representation(2, {[0 -1; 1 0] [1 0; 0 -1]}, true);
dec = D8.naturalRepresentation.irreducible;
U = findUnitary(irrFaith, dec.representation(3));
