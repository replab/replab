function iso1 = niceIsotypicSingleMultiplicity(iso)
% Works on a trivial isotypic component
    if iso.multiplicity ~= 1
        iso1 = replab.DispatchNext('Not a single multiplicity isotypic component');
        return
    end
    irrep = iso.irrep(1);
    irrep1 = replab.nice.niceSubRep(irrep);
    if isequal(irrep1.basis, irrep.basis)
        iso1 = replab.DispatchNext('Recovery did not work');
        return
    end
    iso1 = replab.Isotypic.fromIrreps(iso.parent, {irrep1});
end
