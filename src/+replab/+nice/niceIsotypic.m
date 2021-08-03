function iso1 = niceIsotypic(iso)
    iso1 = [];
    if iso.dimension == 0
        return
    end
    if iso.multiplicity == 1
        irrep1 = replab.nice.niceSubRep(iso.irrep(1));
        if isempty(irrep1)
            return
        end
        iso1 = replab.Isotypic.fromIrreps(iso.parent, {irrep1}, []);
        return
    end
    if iso.trivialDimension == iso.dimension
        sub1 = replab.nice.niceSubRep(iso);
        iso1 = replab.Isotypic.fromTrivialSubRep(sub1);
        return
    end
end
