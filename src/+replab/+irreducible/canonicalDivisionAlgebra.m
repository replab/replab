function res = canonicalDivisionAlgebra(rep, context)
% Computes the representation similar to the given irreducible representation with its division algebra in the RepLAB canonical form
%
% Mutates the argument to add the Frobenius-Schur indicator info.
%
% We return ``res`` such that ``res.parent == rep``.
%
% Args:
%   rep (`+replab.Rep`): Irreducible real representation to find the canonical form of
%
% Returns:
%   `+replab.SimilarRep`: Representation with division algebra in canonical form
    assert(isa(rep, 'replab.Rep'));
    assert(rep.overR);
    assert(isa(context, 'replab.Context'));
    if isempty(rep.frobeniusSchurIndicator)
        rep.frobeniusSchurIndicator = replab.irreducible.frobeniusSchurIndicator(rep, context);
    end
    if ~isequal(rep.isUnitary)
        res = replab.rep.collapse(replab.irreducible.canonicalDivisionAlgebra(rep.unitarize));
        return
    end
    % from then on, we assume unitarity
    switch rep.frobeniusSchurIndicator
      case -1

      case 0
      case 1
      otherwise
        error('Real irreducible representation should have FS indicator -1,0 or 1');
    end
end

end
