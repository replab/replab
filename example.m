addpath(pwd);
addpath([pwd '/external/vpi']);
addpath([pwd '/external/MOxUnit/MOxUnit']);
addpath([pwd '/external/MOxUnit/MOxUnit/util']);
% moxunit_runtests tests -verbose

n = 9;
% construct Sn
Sn = replab.PermutationGroup.fromGenerators([2:n 1], [2 1 3:n]);
g = randperm(n); % or use Sn.sample
disp(sprintf('We factor %s\n', num2str(g)));
w = Sn.factorization(g);
disp(sprintf('to obtain factor %s\n', w.str));
g1 = Sn.evaluateWord(w);
disp(sprintf('evaluated back to %s\n', num2str(g)));
assert(isequal(g, g1));

disp('We use this to construct representations');
rho = Sn.naturalRepresentation('R15');
g = randperm(n);
h = randperm(n);
gh = g(h);
rho.image(g) * rho.image(h)
rho.image(gh)
