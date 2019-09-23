function sub = decomposeUsingCommutant(rep, E1, sampleE, sampleC)
% Decomposes the given representation using eigenvalue decomposition on a generic commutant sample
    tol = replab.Settings.doubleEigTol;
    % Sample a self adjoint commutant element 
    if rep.hasCorrection
        V = diag(rep.D0) * rep.U0;
    else
        V = rep.U0;
    end
    C = full(V * sampleC * V');
    [U D] = replab.rep.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    sub = cell(1, n);
    extra = struct('reducedBlocks', true, 'isIrreducible', true);
    if rep.isExtraFalse('hasTrivialSubspace')
        extra.hasTrivialSubspace = false;
    end
    for i = 1:n
        sub{i} = rep.subRepUnitary(U(:, runs{i})', extra).collapseParent;
    end
end
