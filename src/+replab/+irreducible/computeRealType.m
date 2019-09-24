function t = computeRealType(rep, samples, sub)
% Computes the real type of a real irreducible subrepresentation
%
% Args:
%   rep (replab.Rep): Real representation being decomposed
%   samples (replab.irreducible.OnDemandSamples): Lazy evaluation of various samples for `rep`
%   sub (replab.SubRep): An irreducible subrepresentation of `rep`
%
% Returns:
%   {'R', 'C', 'H'}: Type of the division algebra over the reals
    assert(isequal(rep.field, 'R'));
    assert(sub.isKnownIrreducible);
    assert(sub.parent == rep);
    d = sub.dimension;
    X = sub.U*samples.commutantSample(2)*sub.U'; % can we reuse?
    Xsym = (X+X')/2;
    Xanti = (X-X')/2;
    C = Xsym + 1i * Xanti;
    v = replab.domain.Vectors('C', d).sample;
    v1 = C*v;
    v2 = C*v1;
    tol = replab.Settings.doubleEigTol;
    switch rank([v v1 v2], tol)
      case 1
        t = 'R';
      case 2
        X1 = sub.U*samples.commutantSample(2)*sub.U'; % can we reuse?
        X2 = sub.U*samples.commutantSample(3)*sub.U';
        X3 = sub.U*samples.commutantSample(4)*sub.U';
        H1 = (X1 + X1')/2;
        Hi = (X1 - X1')/2;
        Hj = (X2 - X2')/2;
        Hk = (X3 - X3')/2;
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
