%% Using law checks
%
% This short example demonstrates law checking.
%
% We define the group of permutations acting on 1..10.
n = 10;
eqvFun = @(x, y) isequal(x, y);
sampleFun = @() randperm(n);
composeFun = @(x, y) x(y);
identity = 1:n;
inverseFun = @(x) arrayfun(@(i) find(x == i), 1:10);
S10 = replab.GroupFun('Permutations of 1..10', eqvFun, sampleFun, composeFun, identity, inverseFun)
%%
% We instantiate the group laws, and check them.
S10laws = replab.GroupLaws(S10)
S10laws.check
