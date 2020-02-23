function sub = fullSubRep(rep)
% Creates a full subrepresentation of the given representation with identity basis
%
% Args:
%   rep (`+replab.Rep`): Representaiton
%
% Returns:
%   `+replab.SubRep`: Subrepresentation identical to ``rep``
    d = rep.dimension;
    sub = rep.subRep(speye(d), speye(d));
    assert(isequal(sub.isUnitary, rep.isUnitary));
    sub.trivialDimension = rep.trivialDimension;
    sub.isIrreducible = rep.isIrreducible;
    sub.frobeniusSchurIndicator = rep.frobeniusSchurIndicator;
    sub.isDivisionAlgebraCanonical = rep.isDivisionAlgebraCanonical;
end
