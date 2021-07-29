function f = frobeniusSchurIndicator(rep)
% Computes the Frobenius-Schur indicator of an irreducible subrepresentation
%
% Args:
%   rep (`+replab.Rep`): Representation to compute the Frobenius-Schur indicator of
%
% Returns:
%   integer: Frobenius-Schur indicator
    if rep.overC && rep.inCache('isIrreducible') && rep.isIrreducible
        f = replab.rep.frobeniusSchurIndicator_complex_irreducible(rep);
    elseif isa(rep.group, 'replab.FiniteGroup') && rep.group.order <= replab.globals.autoConjugacyClassesMaximalOrder
        f = replab.rep.frobeniusSchurIndicator_smallFiniteGroup(rep);
    else
        f = replab.rep.frobeniusSchurIndicator_decomposition(rep);
    end
end
