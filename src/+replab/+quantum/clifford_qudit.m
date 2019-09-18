function [G rep] = clifford_qudit(d)
% Returns the Clifford group on a single qudit and its natural representation
%
% G is a permutation group isomorphic to the Clifford group
% rep is the group natural representation on qudits
%
% For now, only d=2 and d=3 are implemented
    switch d
        case 2
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
          rep = G.rep('C', 2, true, {W H S});
      case 3
        % Generators computed using the following in GAP 4
        % omega := E(3);
        % zeta := E(9);
        % H := -E(4)/Sqrt(3)*[[1, 1, 1], [1,omega,omega^2], [1,omega^2,omega]];
        % S := zeta^8 * [[1,0,0],[0,1,0],[0,0,omega]];
        % M := Group(H,S);
        % phi:=IsomorphismPermGroup(M);
        % ListPerm(Image(phi, H));
        % ListPerm(Image(phi, S));
        h = [2 1 5 7 3 10 12 14 16 6 15 20 22 24 11 27 28 30 29 4 34 ...
             13 36 38 40 37 43 45 19 47 48 49 51 21 53 55 54 8 57 25 ...
             59 60 9 63 65 67 69 71 32 61 68 75 39 76 79 73 35 70 84 ...
             83 50 74 44 88 17 82 46 90 18 58 86 94 56 95 52 26 91 89 ...
             23 87 42 64 81 103 93 31 80 66 78 33 96 85 101 72 106 97 ...
             77 107 105 104 92 108 41 98 102 62 100 99];
        s = [3 4 6 8 9 11 13 15 17 18 19 21 23 25 26 2 29 31 32 33 35 ...
             1 37 39 41 42 44 46 38 5 49 50 52 45 54 56 10 36 58 20 60 ...
             61 62 64 66 68 70 72 65 73 74 76 77 78 80 81 82 83 85 75 ...
             86 87 71 14 7 89 43 55 91 92 93 95 22 96 97 98 90 99 84 ...
             100 30 101 102 94 53 16 24 104 79 105 34 28 40 69 12 107 ...
             27 63 48 47 108 67 51 106 88 59 57 103];
        S108 = replab.Permutations(108);
        G = S108.subgroup({h s});
        omega = exp(2i*pi/3);
        zeta = exp(2i*pi/9);
        H = -1i/sqrt(3)*[1 1 1; 1 omega omega^2; 1 omega^2 omega];
        S = zeta^8 * [1 0 0; 0 1 0; 0 0 omega];
        rep = G.rep('C', 3, true, {H S});
      otherwise
        error(sprintf('Dimension d=%d is not supported, only d=2,3', d));
    end
end
