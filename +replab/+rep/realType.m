function t = realType(rep)
% Identifies the type R/C/H of an irreducible real representation
%
%    rep: Real representation
%      t: R for real, C for complex and H for quaternionic
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    C = rep.complexification.commutant.sampleSelfAdjoint;
    v = replab.domain.ComplexVectors(d).sample;
    v1 = C*v;
    v2 = C*v1;
    tol = replab.Settings.doubleEigTol;
    switch rank([v v1 v2], tol)
      case 1
        t = 'R';
      case 2
        H = rep.quaternionification.commutant.sampleSelfAdjoint;
        w = replab.domain.QuaternionVectors(d).sample;
        w1 = H*w;
        w2 = H*w1;
        w3 = H*w2;
        w4 = H*w3;
        switch rank([w w1 w2 w3 w4], tol)
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
