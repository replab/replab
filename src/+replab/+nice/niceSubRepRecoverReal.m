function res = niceSubRepRecoverReal(sub)
% Replaces the complex basis of a subrepresentation by a real basis if possible
    if sub.overR || isreal(sub.basis)
        res = replab.DispatchNext('Basis is already real');
        return
    end
    if ~isequal(sub.isUnitary, true) || ~isequal(sub.parent.isUnitary, true)
        res = replab.DispatchNext('Both the subrepresentation and its parent must be unitary');
        return
    end
    basis = sub.basis;
    d = size(basis, 2);
    v = replab.domain.Vectors('C', d).sample;
    v1 = conj(basis*v);
    if norm(v1 - basis * (basis \ v1)) < replab.globals.doubleEigTol
        % TODO: do we really need that sample? It's there to be sure
        % that generically, the real component is present
        basis1 = orth(real(basis * replab.UnitaryGroup(d).sample));
        sub1 = sub.parent.subRep(basis1, basis1');
        sub1.isIrreducible = sub.isIrreducible;
        sub1.trivialDimension = sub.trivialDimension;
        res = sub1.nice;
    else
        res = replab.DispatchNext('No real basis seems to exist');
        return
    end
end
