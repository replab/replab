function [trivial, nontrivial] = blockTrivialSplit(rep, block, forceNonUnitaryAlgorithms)
% Splits a representation into trivial and nontrivial subrepresentations
%
% It returns two subrepresentations ``trivial`` and ``nontrivial`` whose direct sum is ``rep`` itself.
%
% Args:
%   rep (`+replab.Rep`): Representation to split
%   block (integer(1,\*)): Subset of coordinates describing an existing block of ``rep``
%   forceNonUnitaryAlgorithms (logical): Whether to force the use of a non-unitary algorithm
%
% Returns
% -------
%   trivial: `+replab.SubRep`
%     Trivial subrepresentation, may have ``dimension == 0``
%   nontrivial: `+replab.SubRep`
%     Nontrivial subrepresentation, may have ``dimension == 0``
    D = rep.dimension;
    d = length(block);
    [projTrep, err] = rep.trivialProjector;
    if err >= 1
        error('Precision is not sufficient to compute the trivial dimension.');
    end
    projT = projTrep(block, block);
    projN = eye(d) - projT;
    % dimensions
    dT = round(trace(projT));
    dN = round(trace(projN));
    % find basis
    [IT, PT, pT] = replab.numerical.sRRQR_rank(projT, 1.5, dT);
    [IN, PN, pN] = replab.numerical.sRRQR_rank(projN, 1.5, dN);
    % regularize or apply corrections
    if ~forceNonUnitaryAlgorithms && rep.isUnitary
        PT = IT';
        PN = IN';
    else
        PT(:,pT) = PT;
        PT = (PT*IT)\PT;
        PN(:,pN) = PN;
        PN = (PN*IN)\PN;
    end
    % trivial injection/projection maps
    Itrivial = sparse(D, dT);
    Ptrivial = sparse(dT, D);
    if length(block) > D/2
        Itrivial = full(Itrivial);
        Ptrivial = full(Ptrivial);
    end
    Itrivial(block, :) = IT;
    Ptrivial(:, block) = PT;
    % nontrivial injection/projection maps
    Inontrivial = sparse(D, dN);
    Pnontrivial = sparse(dN, D);
    if length(block) > D/2
        Inontrivial = full(Inontrivial);
        Pnontrivial = full(Pnontrivial);
    end
    Inontrivial(block, :) = IN;
    Pnontrivial(:, block) = PN;
    if ~forceNonUnitaryAlgorithms && rep.isUnitary
        nontrivial = rep.subRep(Inontrivial, 'projection', Pnontrivial, 'isUnitary', true, 'trivialDimension', 0);
        trivial = rep.subRep(Itrivial, 'projection', Ptrivial, 'isUnitary', true, 'trivialDimension', dT);
    else
        nontrivial = rep.subRep(Inontrivial, 'projection', Pnontrivial, 'trivialDimension', 0);
        trivial = rep.subRep(Itrivial, 'projection', Ptrivial, 'isUnitary', true, 'trivialDimension', dT);
    end
end
