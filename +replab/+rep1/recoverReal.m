function U1 = recoverReal(sub)
% Replaces the basis sub.U by a real basis if possible
    assert(isequal(sub.field, 'C'));
    U = sub.U;
    d = size(U, 2);
    v = replab.domain.ComplexVectors(d).sample;
    v1 = conj(U*v);
    if norm(v1 - U * (U \ v1)) < replab.Settings.doubleEigTol
        U1 = orth(real(U * replab.domain.UnitaryMatrices(d)));
    else
        U1 = [];
    end
end
