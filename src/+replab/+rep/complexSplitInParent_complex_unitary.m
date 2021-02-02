function irreps = complexSplitInParent_complex_unitary(sub, iterator)
% Decomposes fully a complex subrepresentation into irreducible subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Unitary subrepresentation to split further
%   iterator (`+replab.+domain.SamplesIterator`): Iterator in the sequence of parent commutant samples
%x
% Returns:
%   cell(1,\*) of `.SubRep`: Irreducible subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.isUnitary);
    assert(sub.overC);
    tol = replab.globals.doubleEigTol;
    d = sub.dimension;
    % force a dense matrix to have eig behave well
    S = full(sub.projection('double/sparse') * iterator.next * sub.injection('double/sparse'));
    X = (S + S')/2;
    [U, D] = eig(X);
    D = reshape(diag(D), [1 d]);
    D = D(:).';
    P = replab.Partition.fromApproximateVector(D, tol);
    n = P.nBlocks;
    if n == 1
        sub.cache('isIrreducible', true, '==');
        irreps = {sub};
    else
        irreps = cell(1, n);
        for i = 1:n
            basis = U(:, P.block(i));
            I = sub.injection('double/sparse') * basis;
            P = basis' * sub.projection('double/sparse');
            irreps{i} = sub.parent.subRep(I, 'projection', P, 'isUnitary', true, 'isIrreducible', true);
        end
    end
    if sub.inCache('trivialDimension') && sub.trivialDimension == 0
        for i = 1:length(irreps)
            irreps{i}.cache('trivialDimension', 0, '==');
        end
    end
end
