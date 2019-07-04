function [Utrivial Urest] = extractTrivial(rep)
    d = rep.dimension;
    T = replab.trivial(rep);
    O = replab.rep.orbits(rep);
    U = zeros(d, 0);
    D = zeros(0, 1);
    for i = 1:O.nBlocks
        v = zeros(d, 1);
        b = O.block(i);
        v(b) = 1;
        if T.isInvariant(v)
            v = v;
            U = [U v];
            D(1, i) = 1/length(b);
        end
    end
    nRest = d - size(U, 2);
    switch rep.field
      case 'R'
        M = replab.domain.RealMatrices(d, nRest);
      case 'C'
        M = replab.domain.ComplexMatrices(d, nRest);
      case 'R'
        M = replab.domain.QuaternionMatrices(d, nRest);
    end
    S = M.sample;
    if size(U, 2) > 0
        S = S - U*diag(D)*U'*S;
    end
    S = T.project(S);
    if M.eqv(S, zeros(d, nRest))
        Utrivial = U;
    else
        Utrivial = [U orth(S)];
    end
    if size(Utrivial, 2) > 0
        Urest = null(Utrivial');
    else
        Urest = speye(d);
    end
    assert(size(Utrivial, 1) == d);
    assert(size(Urest, 1) == d);
    assert(size(Utrivial, 2) + size(Urest, 2) == d);
end
