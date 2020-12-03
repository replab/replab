replab.globals.verbosity(2); % try 1 or 2
n = 5;
G = replab.S(n);
% get the natural representation
nat = G.naturalRep;
triv = nat.subRep(ones(n, 1));
% split the trivial representation and the standard representation using Maschke's theorem
std = nat.maschke(triv);
% construct a tensor power
pow4 = nat.tensorPower(4);
% and a subrepresentation of that tensor power
sub = pow4.subRep(kron(kron(triv.injection, std.injection), kron(triv.injection, std.injection)));
% construct a representation with tons of noise
subNoisy = sub.withNoise(0.01);
% the error bound is horrific
fprintf('Error bound with noise: %e\n', subNoisy.errorBound);
% but we can recover a pretty precise subrepresentation
tic; subRefined = subNoisy.refine; t = toc;
fprintf('Time elapsed: %f\n', t);
fprintf('Error bound after refinement: %e\n', subRefined.errorBound);
fprintf('Now let us use the large scale variant\n');
% but we can recover a pretty precise subrepresentation
tic; subRefined = subNoisy.refine('largeScale', true, 'nSamples', 5); t = toc;
fprintf('Time elapsed: %f\n', t);
fprintf('Error bound after refinement (large scale): %e\n', subRefined.errorBound);
% the large scale algorithm is slower, but it is the one that works well for super large
% dimensions (d > 5000)