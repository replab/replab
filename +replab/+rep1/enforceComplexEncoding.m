function [T V] = enforceComplexEncoding(rep)
    assert(isequal(rep.field, 'R'));
    d = rep.dimension;
    A = rep.complexification.commutant.sampleSelfAdjoint;
    v1 = replab.domain.ComplexVectors(d).sample;
    v1 = v1/norm(v1);
    w1p = A*v1;
    alpha1 = real(w1p'*v1);
    w1 = w1p - alpha1*v1;
    beta2 = norm(w1);
    v2 = w1/beta2;
    w2p = A*v2;
    alpha2 = real(w2p'*v2);
    w2 = w2p - alpha2*v2 - beta2*v1;
    norm(w2)
    T = [alpha1 beta2
         beta2 alpha2];
    V = [v1 v2];
end
