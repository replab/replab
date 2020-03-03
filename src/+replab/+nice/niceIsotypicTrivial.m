function iso1 = niceIsotypicTrivial(iso)
% Works on a trivial isotypic component
    if ~isequal(iso.trivialDimension, iso.dimension)
        iso1 = replab.DispatchNext('Not a trivial isotypic component');
        return
    end
    sub1 = replab.nice.niceSubRep(iso);
    if isequal(sub1.basis, iso.basis)
        iso1 = replab.DispatchNext('Recovery did not work');
        return
    end
    iso1 = replab.Isotypic.fromTrivialSubRep(iso, sub1);
end
