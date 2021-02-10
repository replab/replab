function subs = absoluteSplitInParent_complex_unitary(sub, sample)
% Decomposes a complex unitary subrepresentation into subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Unitary subrepresentation to split further
%   sample (double(\*,\*)): Sample of ``sub.parent.commutant``
%x
% Returns:
%   cell(1,\*) of `.SubRep`: Subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.isUnitary);
    assert(sub.overC);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    subI = sub.injection('double/sparse');
    subP = sub.projection('double/sparse');
    S = full(subP * sample * subI);
    X = (S + S')/2;
    [U, D] = eig(X);
    D = reshape(diag(D), [1 d]);
    P = replab.Partition.fromApproximateVector(D, tol);
    blocks = P.blocks;
    n = P.nBlocks;
    if n == 1
        sub.cache('isIrreducible', true, '==');
        subs = {sub};
    else
        subs = cell(1, n);
        for i = 1:n
            basis = U(:, blocks{i});
            I = subI * basis;
            if all(subI == subP')
                P = I';
            else
                P = basis' * subP;
            end
            subs{i} = sub.parent.subRep(I, 'projection', P, 'isUnitary', true);
        end
    end
end
