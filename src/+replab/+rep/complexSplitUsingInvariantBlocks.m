function irreps = complexSplitUsingInvariantBlocks(rep)
% Splits a representation into irreducible subrepresentations using existing block structure
%
% First extracts the trivial subrepresentations before splitting the nontrivial invariant subspace;
% this is done block by block.
%
% Args:
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Irreducible subrepresentations of ``rep``
    partition = rep.invariantBlocks;
    t = cputime;
    D = rep.dimension;
    P = speye(D);
    [P1, E1] = rep.trivialRowSpace.project(P, 'double');
    [P2, E2] = rep.trivialColSpace.project(P1, 'double');
    replab.msg(2, 'Time (trivial subspace projection): %2.2f s', cputime - t);
    if E1 + E2 >= 1
        error('Representation is not precise enough to compute the trivial dimension.');
    end
    Itriv = sparse(D, 0);
    Ptriv = sparse(0, D);
    nontrivialIrreps = cell(1, 0);
    for i = 1:partition.nBlocks
        blk = partition.block(i);
        d = length(blk);
        % projectors on trivial and nontrivial space
        projT = P2(blk, blk);
        projN = eye(d) - projT;
        % dimensions
        dT = round(trace(projT));
        dN = round(trace(projN));
        [IT, PT, pT] = replab.numerical.sRRQR_rank(projT, 1.5, dT);
        [IN, PN, pN] = replab.numerical.sRRQR_rank(projN, 1.5, dN);
        % regularize or apply corrections
        if rep.knownUnitary
            PT = IT';
            PN = IN';
        else
            PT(:,pT) = PT;
            PT = (PT*IT)\PT;
            PN(:,pN) = PN;
            PN = (PN*IN)\PN;
        end
        ITnew = sparse(D, dT);
        PTnew = sparse(dT, D);
        ITnew(blk, :) = IT;
        PTnew(:, blk) = PT;
        % update basis of trivial space
        Itriv = [Itriv ITnew];
        Ptriv = [Ptriv; PTnew];
        INnew = sparse(D, dN);
        PNnew = sparse(dN, D);

        INnew(blk, :) = IN;
        PNnew(:, blk) = PN;
        if rep.knownUnitary
            subN = rep.subRep(INnew, 'projection', PNnew, 'isUnitary', true, 'trivialDimension', 0);
        else
            subN = rep.subRep(INnew, 'projection', PNnew, 'trivialDimension', 0);
        end
        nontrivialIrreps = horzcat(nontrivialIrreps, subN.complexSplitInParent);
    end
    if rep.knownUnitary
        subT = rep.subRep(Itriv, 'projection', Ptriv, 'isUnitary', true);
    else
        subT = rep.subRep(Itriv, 'projection', Ptriv);
    end
    trivialIrreps = replab.Isotypic.fromTrivialSubRep(subT).irreps;
    irreps = horzcat(trivialIrreps, nontrivialIrreps);
end
