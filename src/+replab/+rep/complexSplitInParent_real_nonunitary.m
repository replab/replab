function irreps = complexSplitInParent_real_nonunitary(sub, iterator)
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
    error('stop')
    % TODO: implement approximate vector with complex numbers
    %P = replab.Partition.fromApproximateVector(real(D), tol);
    %blocks = P.blocks;
    %n = P.nBlocks;
    %if n == 1
    %    sub.cache('isIrreducible', true, '==');
    %    irreps = {sub};
    %else
    %    irreps = cell(1, n);
    %    for i = 1:n
    %        inj1 = U(:, blocks{i});
    %        prj1 = V(blocks{i}, :);
    %        prj1 = (prj1*inj1)\prj1;
    %        I = sub.injection('double/sparse') * inj1;
    %        P = prj1 * sub.projection('double/sparse');
    %        irreps{i} = sub.parent.subRep(I, 'projection', P, 'isIrreducible', true);
    %    end
    %end
    %if sub.inCache('trivialDimension') && sub.trivialDimension == 0
    %    for i = 1:length(irreps)
    %        irreps{i}.cache('trivialDimension', 0, '==');
    %    end
    %end
end
