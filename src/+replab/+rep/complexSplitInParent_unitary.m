function irreps = complexSplitInParent_unitary(sub, iterator)
% Decomposes fully a subrepresentation into irreducible subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Unitary subrepresentation to split further
%   iterator (`+replab.+domain.SamplesIterator`): Iterator in the sequence of parent commutant samples
%x
% Returns:
%   cell(1,\*) of `.SubRep`: Irreducible subrepresentations with their ``.parent`` set to the ``.parent`` of ``sub``
    assert(sub.isUnitary);


    tol = replab.Parameters.doubleEigTol;
    % extract nontrivial representations by sampling the commutant
    C = full(self.projection('double/sparse') * self.parent.commutant.sample('double') * self.injection('double/sparse'));
    % the eigenspace decomposition is the basis of the numerical decomposition
    % V'*C*V = D
    [U1 D] = replab.numerical.sortedEig((C + C')/2, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    % find repeated eigenvalues
    runs = replab.UndirectedGraph.fromAdjacencyMatrix(mask).connectedComponents.blocks;
    n = length(runs);
    if n == 1
        self.cache('isIrreducible', true, '==');
        irreps = {self};
    else
        irreps = cell(1, n);
        for i = 1:n
            basis = U1(:, runs{i});
            I = self.injection('double/sparse') * basis;
            P = basis' * self.projection('double/sparse');
            irreps{i} = self.parent.subRep(I, 'projection', P, 'isUnitary', true, 'isIrreducible', true);
        end
    end
    if self.inCache('trivialDimension') && self.trivialDimension == 0
        for i = 1:length(irreps)
            irreps{i}.cache('trivialDimension', 0, '==');
        end
    end
end
