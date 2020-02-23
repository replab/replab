function res = collapseSubRepSubRep(rep)
% Given a subrepresentation of a subrepresentation, returns a subrepresentation
%
% Args:
%   rep (`+replab.SubRep`): A subrepresentation with ``rep.parent`` an instance of `+replab.SubRep` too
%
% Returns:
%   `+replab.SubRep`: A subrepresentation of ``rep.parent.parent``
    assert(isa(rep, 'replab.SubRep'));
    parent = rep.parent;
    assert(isa(parent, 'replab.SubRep'));
    newB_internal = parent.B_internal * rep.B_internal;
    newE_internal = rep.E_internal * parent.E_internal;
    res = parent.parent.subRep(newB_internal, newE_internal);
    res.trivialDimension = rep.trivialDimension;
    res.isIrreducible = rep.isIrreducible;
    res.frobeniusSchurIndicator = rep.frobeniusSchurIndicator;
    res.isDivisionAlgebraCanonical = rep.isDivisionAlgebraCanonical;
end
