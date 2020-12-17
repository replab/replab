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
%   context (`+replab.Context`): Sampling context
%
% Returns:
%   `+replab.SimilarRep`: Representation with division algebra in canonical form
    assert(isa(rep, 'replab.Rep') && rep.isIrreducible);
    assert(isa(context, 'replab.Context'));

    if rep.overC || rep.frobeniusSchurIndicator == 1
        % if over complex field, or real-type representation, no need to change
        res = replab.SimilarRep.identical(rep);
        return
    end

    % from then on, we assume unitarity
    switch rep.frobeniusSchurIndicator
      case -2
        % quaternion
        res = replab.irreducible.enforceQuaternionEncoding(rep, context);
      case 0
        % complex
        res = replab.irreducible.enforceComplexEncoding(rep, context);
      otherwise
        error('Real irreducible representation should have FS indicator -2, 0 or 1');
    end
end
