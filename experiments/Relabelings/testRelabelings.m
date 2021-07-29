run ../../replab_init.m
% We construct the relabelings for the scenario with nParties = 2, mInputs = 2, kOutputs = 2
n = 2;
m = 2;
k = 2;
mu = relabelings(n, m, k, 'Pabxy');
% mu is an isomorphism with a structured source group, and an unstructured target group
mu.source % this is a double wreath product group
mu.target % this is a permutation group representing the action on P(a,b,x,y)

% so basically, you can perform everything using mu.target, and discard mu.source and mu itself

% We construct the relabeling action on P(a,b,x,y) by hand to verify the construction
Iabxy = reshape(1:16, [2 2 2 2]);
IPermParties = reshape(permute(Iabxy, [2 1 4 3]), 1, []);
IPermAliceInput = Iabxy(:,:,[2 1],:);
IPermAliceOutputForFirstInput = Iabxy;
IPermAliceOutputForFirstInput(:,:,1,:) = IPermAliceOutputForFirstInput([2 1],:,1,:);

% The corresponding elements of the double wreath product
wPermParties = {[2, 1], { ... % Permutation of parties
    {[1, 2], ... % Permutation of Alice inputs
     {[1, 2], [1, 2]}}, ... % Permutation of Alice outputs for x=1 and x=2
    {[1, 2], ... % Permutation of Bob inputs
     {[1, 2], [1, 2]}} ... % Permutation of Bob outputs for y=1 and y=2
                   }};
wPermAliceInput = {[1, 2], { ...
    {[2, 1], {[1, 2], [1, 2]}}, ...
    {[1, 2], {[1, 2], [1, 2]}}
                   }};
wPermAliceOutputForFirstInput = {[1, 2], { ...
    {[1, 2], {[2, 1], [1, 2]}}, ...
    {[1, 2], {[1, 2], [1, 2]}}
                   }};
% Verification
assert(all(mu.imageElement(wPermParties) == reshape(IPermParties, 1, [])));
assert(all(mu.imageElement(wPermAliceInput) == reshape(IPermAliceInput, 1, [])));
assert(all(mu.imageElement(wPermAliceOutputForFirstInput) == reshape(IPermAliceOutputForFirstInput, 1, [])));

% Construct the CHSH inequality
I = zeros(2,2,2,2);
for a = 0:1
    for b = 0:1
        for x = 0:1
            for y = 0:1
                I(a+1,b+1,x+1,y+1) = (-1)^(a+b+x*y);
            end
        end
    end
end
% Construct the relabeling subgroup that leaves CHSH invariant
G = mu.target; % group of all relabelings as permutations
Iflat = reshape(I, 1, []); % argument of permutation group methods must be row vectors
GCHSH = G.vectorStabilizer(Iflat);
WCHSH = mu.preimageGroup(GCHSH) % for fun, we can obtain the structured subgroup to get insights

% now, let us find the canonical representative of CHSH
Ican = G.vectorFindLexMinimal(Iflat, GCHSH); % we pass GCHSH (optional 2nd arg) so we avoid recomputation

% now, let us construct all the relabelings of CHSH
% we know that g(I) = I for all g in GCHSH but this is not the case for all g in G
%
% thus we want to split the elements of G in equivalence classes
% g1 ~ g2 if there exists h such that g1 = g2 h
%
% the correct way is to enumerate the left cosets of GCHSH in G
cosets = G.leftCosetsOf(GCHSH);
assert(cosets.nElements == 8); % number of representatives
explicitCosets = cosets.elements;
reps = zeros(length(explicitCosets), length(Iflat));
for i = 1:length(explicitCosets)
    g = explicitCosets{i}.representative;
    gInv = G.inverse(g);
    reps(i,:) = Iflat(gInv); % the action of a permutation g on a vector I is given by I(inv(g))
end