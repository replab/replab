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
rho = Sn.naturalRepresentation;
g = randperm(n);
h = randperm(n);
gh = g(h);
rho.image(g) * rho.image(h)
rho.image(gh)

rho.nIsotypicComponents
[subrho1 U1 dim1 mul1] = rho.isotypicComponent(1)

rho.nIsotypicComponents
[subrho2 U2 dim2 mul2] = rho.isotypicComponent(2)
