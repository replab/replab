function t = realType(rep)
% Identifies the type R/C/H of an irreducible real representation
%
%    rep: Real representation
%      t: R for real, C for complex and H for quaternionic
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    C = rep.complexification.commutant.sampleSelfAdjoint;
    v = replab.domain.Vectors('C', d).sample;
    v1 = C*v;
    v2 = C*v1;
    tol = replab.Parameters.doubleEigTol;
    switch rank([v v1 v2], tol)
      case 1
        t = 'R';
      case 2
        H1 = rep.commutant.sampleSelfAdjoint;
        Hi = rep.commutant.sampleSkewAdjoint;
        Hj = rep.commutant.sampleSkewAdjoint;
        Hk = rep.commutant.sampleSkewAdjoint;
        H = replab.quaternion.H(H1, Hi, Hj, Hk);
        w = replab.quaternion.Vectors(d).sample;
        w1 = H*w;
        w2 = H*w1;
        w3 = H*w2;
        w4 = H*w3;
        W = [w.r w1.r w2.r w3.r w4.r
             w.i w1.i w2.i w3.i w4.i
             w.j w1.j w2.j w3.j w4.j
             w.k w1.k w2.k w3.k w4.k];
        switch rank(W, tol)
          case 2
            t = 'C';
          case 4
            t = 'H';
          otherwise
            error('Problem when identifying C/H irrep type'); 
        end
      otherwise
        error('Problem when identifying R/C/H irrep type'); 
    end
end
