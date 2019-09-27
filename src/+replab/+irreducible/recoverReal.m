function U1 = recoverReal(U)
% Replaces the basis U by a real basis if possible
    U = U';
    d = size(U, 2);
    v = replab.domain.Vectors('C', d).sample;
    v1 = conj(U*v);
    if norm(v1 - U * (U \ v1)) < replab.Settings.doubleEigTol
        U1 = orth(real(U * replab.UnitaryGroup(d).sample));
        U1 = U1'; % return to row vector convention
    else
        U1 = [];
    end
end
