%% Example of replab usage
%
% Plottons la fonction
% $\sin(x)$
%

x =[0:0.01:1];
plot(x,sin(x));

% Run first replab_addpaths;

n = 4;
% construct Sn
Sn = replab.PermutationGroup.fromGenerators({[2:n 1] [2 1 3:n]});
% sample permutation; alternative is to use Sn.sample or Sn.sampleUniformly
g = randperm(n);
disp(sprintf('We factor %s\n', num2str(g)));
w = Sn.factorization(g);
disp(sprintf('to obtain factor %s\n', w.str));
g1 = Sn.evaluateWord(w);
disp(sprintf('evaluated back to %s\n', num2str(g)));
assert(isequal(g, g1));

% ceci est un commentaire

%%
% Un peu de texte ici

%% Nouvelle section
disp('We use this to construct representations and evaluate them');
rho = Sn.naturalRepresentation;
g = randperm(n);
h = randperm(n);
gh = g(h);
rho.image(g) * rho.image(h)
rho.image(gh)

%%
% Un peu de texte par là

[subrho1 U1] = rho.isotypic.component(1)

[subrho2 U2] = rho.isotypic.component(2)

M = rho.centralizerAlgebra.project(rand(n,n))
blocks = rho.isotypic.blocksOf(M)
blocks{1}
blocks{2}
