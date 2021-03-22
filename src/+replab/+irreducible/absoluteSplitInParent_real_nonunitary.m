function subs = absoluteSplitInParent_real_nonunitary(sub, sample)
% Decomposes a real subrepresentation into subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to split further
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%x
% Returns:
%   cell(1,\*) of `.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.overR);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    S = full(sub.projection('double/sparse') * sample * sub.injection('double/sparse'));
    [U, D, V] = replab.numerical.realeig(S);
    V = V';
    Dreal = reshape(diag(D), 1, []);
    P = replab.Partition.fromApproximateVector(real(Dreal), tol);
    blocks = P.blocks;
    n = P.nBlocks;
    subs = cell(1, n);
    for i = 1:n
        blk = blocks{i};
        inj1 = U(:, blk);
        prj1 = V(blk, :);
        prj1 = (prj1*inj1)\prj1;
        I = sub.injection('double/sparse') * inj1;
        P = prj1 * sub.projection('double/sparse');
        isReal = all(abs(diag(D(blk(1:end-1),blk(2:end)))) <= tol);
        if isReal
            subs{i} = sub.parent.subRep(I, 'projection', P);
        else
            subs{i} = sub.parent.subRep(I, 'projection', P, 'divisionAlgebraName', 'C->R');
        end
    end
end
