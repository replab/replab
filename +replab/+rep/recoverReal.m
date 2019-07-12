function U1 = recoverReal(sub)
% Replaces the basis U by a real basis if possible
    U = sub.U'; % use column vectors
    d = size(U, 2);
    v = replab.domain.ComplexVectors(d).sample;
    v1 = conj(U*v);
    if norm(v1 - U * (U \ v1)) < replab.Settings.doubleEigTol
        U1 = orth(real(U * replab.domain.UnitaryMatrices(d)));
        U1 = U1'; % return to row vector convention
    else
        U1 = [];
    end
end
