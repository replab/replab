% The GHZ state |GHZ> = |000> + |111> is invariant (U |GHZ> = |GHZ>) under a family of matrices U
%
%      [a0  0]      [b0  0]     [c0  0]
% U =  [ 0 a1]  (x) [ 0 b1] (x) [ 0 c1]
%
% who act on the state space ``C^8``. The eight coefficients, in order, are those of
% ``|000>``, ``|001>``, ``|010>``, ``|011>``, ``|100>``, ``|101>``, ``|110>``, ``|111>``.
%
% For the invariance to hold, the unit complex numbers a0,a1,b0,b1,c0,c1 obey the equation a0*b0*c0 = a1*b1*c1 = 1.
%
% The state is also invariant under permutation of subsystems (a group of order 3! = 6), and permutation of the levels
% (a group of order 2).

% Let us define the continuous connected group.
%
% T6 is the torus group with 6 elements, that we name according to our scenario:
T6 = replab.T(6).withNames({'a0' 'b0' 'c0' 'a1' 'b1' 'c1'});

% We define the subgroup obeying the required equations
T = T6.subgroupWith('a0*b0*c0 = 1', 'a1*b1*c1 = 1');

% Now, how does that group act on the state space ``C^8``?
%
% We construct the representation whose image is the matrix U above:
Trep_eqs = {};
Trep = T.diagonalRepWith('a0 b0 c0', ...
                         'a0 b0 c1', ...
                         'a0 b1 c0', ...
                         'a0 b1 c1', ...
                         'a1 b0 c0', ...
                         'a1 b0 c1', ...
                         'a1 b1 c0', ...
                         'a1 b1 c1');

% We construct now the discrete part, by writing how the subsystem and level permutations affect
% the elements of the continuous connected part.

% Permutation of AB
gAB = [2 1 3 5 4 6];
% Permutation of AC
gAC = [3 2 1 6 5 4];
% Permutation of BC
gBC = [1 3 2 4 6 5];

% Permutation of the two levels
gL = [4 5 6 1 2 3];

% The finite, discrete group part
F = replab.PermutationGroup.of(gAB, gAC, gBC, gL);

% How does this act on the state space?
% Remember 1,2,3,4,5,6,7,8 enumerates
% ``|000>``, ``|001>``, ``|010>``, ``|011>``, ``|100>``, ``|101>``, ``|110>``, ``|111>``.

% to simplify things, we first write a morphism from F to S(8), with the following images
hAB = [1 2 5 6 3 4 7 8];
hAC = [1 5 3 7 2 6 4 8];
hBC = [1 3 2 4 5 7 6 8];
hL = [8 7 6 5 4 3 2 1];

% here is the morphism
mu = F.morphismByImages(replab.S(8), 'preimages', {gAB, gAC, gBC, gL}, 'images', {hAB, hAC, hBC, hL});

% and we construct a representation by composing the morphism with the permutation matrix representation of S(8)
Frep = mu.andThen(replab.S(8).naturalRep);

G = T.semidirectProductFromRep(F.naturalRep);
rep = G.semidirectProductRep(Frep.complexification, Trep);
