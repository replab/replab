function irreps = absoluteSplitInParent_real_nonunitary(sub, iterator)
% Decomposes fully a real subrepresentation into irreducible subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to split further
%   iterator (`+replab.+domain.SamplesIterator`): Iterator in the sequence of parent commutant samples
%x
% Returns:
%   cell(1,\*) of `.SubRep`: Irreducible subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.overR);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    S = full(sub.projection('double/sparse') * iterator.next * sub.injection('double/sparse'));
    [U, D, V] = replab.numerical.realeig(S);
    V = V';
    Dreal = reshape(diag(D), 1, []);
    P = replab.Partition.fromApproximateVector(real(Dreal), tol);
    blocks = P.blocks;
    n = P.nBlocks;
    irreps = cell(1, n);
    for i = 1:n
        blk = blocks{i};
        inj1 = U(:, blk);
        prj1 = V(blk, :);
        prj1 = (prj1*inj1)\prj1;
        I = sub.injection('double/sparse') * inj1;
        P = prj1 * sub.projection('double/sparse');
        isReal = all(abs(diag(D(blk(1:end-1),blk(2:end)))) <= tol);
        if isReal
            irreps{i} = sub.parent.subRep(I, 'projection', P, 'isIrreducible', true, 'frobeniusSchurIndicator', 1);
        else
            irreps{i} = sub.parent.subRep(I, 'projection', P, 'divisionAlgebraName', 'complex');
        end
    end
    if sub.inCache('trivialDimension') && sub.trivialDimension == 0
        for i = 1:length(irreps)
            irreps{i}.cache('trivialDimension', 0, '==');
        end
    end
end
