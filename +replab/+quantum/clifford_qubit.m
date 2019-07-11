function [G rep] = clifford_qubit
% Returns the Clifford group on a single qubit and its natural representation
%
% G is a permutation group isomorphic to the Clifford group
% rep is the group natural representation on qubits
    S48 = replab.Permutations(48);
    w = [2 5 1 6 11 12 3 4 13 14 20 21 22 23 7 8 9 10 24 15 30 ...
         31 32 33 16 17 18 19 34 25 41 42 43 39 26 27 28 29 44 38 ...
         35 36 48 47 37 40 46 45];
    h = [4 6 8 1 12 2 16 3 14 17 21 5 23 9 25 7 10 26 29 30 11 ...
         32 13 34 15 18 35 38 19 20 42 22 39 24 27 41 40 28 33 37 ...
         36 31 44 43 46 45 48 47];
    s = [1 2 3 9 5 13 7 17 19 4 11 22 24 6 15 26 28 8 10 20 31 ...
         33 12 14 35 37 16 18 39 41 43 21 23 44 45 25 27 34 47 29 ...
         48 30 32 46 36 38 40 42];
    G = S48.subgroup({w h s});
    E8 = exp(1i*pi/4);
    E4 = 1i;
    W = [E8 0; 0 E8];
    H = [1 1; 1 -1]/sqrt(2);
    S = [1 0; 0 E4];
    rep = G.rep('C', 2, {W H S});
end
