function res = niceSubRepRecoverInteger(sub)
% Replaces the real basis of a subrepresentation by an integer basis if possible
    if ~isreal(sub.basis)
        res = replab.DispatchNext('Basis needs to be real');
        return
    end
    if isequal(sub.parent.isUnitary, true)
        intBasisOpt = replab.nice.integerRowSpan(sub.basis');
        if isempty(intBasisOpt)
            res = replab.DispatchNext('No integer basis found');
            return
        end
        res = sub.parent.subRep(intBasisOpt');
        res.trivialDimension = sub.trivialDimension;
        res.isIrreducible = sub.isIrreducible;
    else
        res = replab.DispatchNext('Parent representation must be unitary');
        return
    end
end
