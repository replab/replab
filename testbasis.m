M2 = blkdiag([0 -1; -1 0], eye(4));
M1 = [[-1 1 -1 1 -1] 1
      -eye(5) zeros(5,1)];
rep2 = G.repByImages('R', 6, {M1 M2}, {inv(M1) inv(M2)});
rep1 = G.definingRep;