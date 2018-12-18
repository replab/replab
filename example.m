% Run first replab_addpaths;

n = 9;
% construct Sn
Sn = replab.PermutationGroup.fromGenerators([2:n 1], [2 1 3:n]);
% sample permutation; alternative is to use Sn.sample or Sn.sampleUniformlyx
g = randperm(n);
disp(sprintf('We factor %s\n', num2str(g)));
w = Sn.factorization(g);
disp(sprintf('to obtain factor %s\n', w.str));
g1 = Sn.evaluateWord(w);
disp(sprintf('evaluated back to %s\n', num2str(g)));
assert(isequal(g, g1));

disp('We use this to construct representations and evaluate them');
rho = Sn.naturalRepresentation('R15');
g = randperm(n);
h = randperm(n);
gh = g(h);
rho.image(g) * rho.image(h)
rho.image(gh)
disp('Random group algebra element')
rho.sampleGroupAlgebra