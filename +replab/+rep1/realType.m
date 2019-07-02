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
    switch rank([v v1 v2])
      case 1
        t = 'R';
      case 2
        H = rep.quaternionification.commutant.sampleSelfAdjoint;
        v = replab.domain.QuaternionVectors(d).sample;
        v1 = H*v;
        v2 = H*v1;
        v3 = H*v2;
        v4 = H*v3;
        switch rank([v v1 v2 v3 v4])
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
