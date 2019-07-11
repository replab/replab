function sub = decomposeUsingCommutant(rep)
% Decomposes the given representation into irreducible representations
%
% Uses eigenvalue decomposition on a generic commutant sample
    assert(isa(rep, 'replab.Rep'));
    tol = replab.Settings.doubleEigTol;
    % sample a Hermitian 
    C = full(rep.commutant.sampleSelfAdjoint);
    [U D] = replab.rep.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    sub = cell(1, n);
    for i = 1:n
        sub{i} = rep.subRep(U(:, runs{i})');
    end
end
