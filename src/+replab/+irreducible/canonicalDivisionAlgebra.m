function res = canonicalDivisionAlgebra(rep, context)
% Computes the representation similar to the given irreducible representation with its division algebra in the RepLAB canonical form
%
% Mutates the argument to add the Frobenius-Schur indicator info.
%
% We return ``res`` such that ``res.parent == rep``.
%
% If ``rep`` is over the complex field, or the real representation type is "real", then we return
% a similar representation with an identity change of basis.
%
% Args:
%   rep (`+replab.Rep`): Irreducible representation to find the canonical form of
%
% Returns:
%   `+replab.SimilarRep`: Representation with division algebra in canonical form
    assert(isa(rep, 'replab.Rep'));
    assert(isa(context, 'replab.Context'));
    if rep.overC
        res = replab.SimilarRep.identical(rep);
        return
    end
    if isempty(rep.frobeniusSchurIndicator)
        rep.frobeniusSchurIndicator = replab.irreducible.frobeniusSchurIndicator(rep, context);
    end
    if ~isequal(rep.isUnitary, true)
        res = replab.rep.collapse(replab.irreducible.canonicalDivisionAlgebra(rep.unitarize, context));
        return
    end
    % from then on, we assume unitarity
    switch rep.frobeniusSchurIndicator
      case -1
        % quaternion
        res = replab.irreducible.enforceQuaternionEncoding(rep);
      case 0
        % complex
        res = replab.irreducible.enforceComplexEncoding(rep);
      case 1
        % real
        res = replab.SimilarRep.identical(rep);
      otherwise
        error('Real irreducible representation should have FS indicator -1,0 or 1');
    end
end
