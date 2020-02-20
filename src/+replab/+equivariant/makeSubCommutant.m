function e = makeSubCommutant(repC, repR, special)
    if ~isequal(special, 'commutant')
        e = replab.DispatchNext('Must be a commutant space');
        return
    end
    if ~isa(repC, 'replab.SubRep')
        e = replab.DispatchNext('Must be a subrepresentation');
        return
    end
    e = repC.parent.commutant.subEquivariant(repC, repR, special);
end
