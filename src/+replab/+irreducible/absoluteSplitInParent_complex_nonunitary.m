function subs = absoluteSplitInParent_complex_nonunitary(sub, sample)
% Decomposes a complex subrepresentation into subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to split further
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.overC);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    S = full(sub.projection('double/sparse') * sample * sub.injection('double/sparse'));
    [U, D, V] = eig(S);
    V = V';
    D = reshape(diag(D), [1 d]);
    % TODO: implement approximate vector with complex numbers
    P = replab.Partition.fromApproximateVector(real(D), tol);
    blocks = P.blocks;
    n = P.nBlocks;
    if n == 1
        subs = {sub};
    else
        subs = cell(1, n);
        for i = 1:n
            blk = blocks{i};
            inj1 = U(:, blk);
            prj1 = V(blk, :);
            prj1 = (prj1*inj1)\prj1;
            I = sub.injection('double/sparse') * inj1;
            P = prj1 * sub.projection('double/sparse');
            subs{i} = sub.parent.subRep(I, 'projection', P);
        end
    end
end
