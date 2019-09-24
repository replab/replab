function sub = splitUsingCommutant(rep, samples, U)
    d = rep.dimension;
    if nargin < 3
        U = speye(d);
    end
    dU = size(U, 1);
    tol = replab.Settings.doubleEigTol;
    if rep.overR
        trivialIrrepInfo = replab.IrrepInfo('1', 'R', []);
    else
        trivialIrrepInfo = replab.IrrepInfo('1', [], []);
    end
    % extract trivial representations
    S = U*samples.trivialSample(1)*U';
    trivials = {};
    if ~replab.isNonZeroMatrix(S, replab.Settings.doubleEigTol) % magic epsilon
        Utrivial = orth(S)';
        Urest = null(Utrivial)';
        assert(size(Utrivial, 2) == dU);
        assert(size(Urest, 2) == dU);
        assert(size(Urest, 1) + size(Utrivial, 1) == dU);
        for i = 1:size(Utrivial, 1)
            trivials{1,i} = rep.subRepUnitary(Utrivial * U, [], trivialIrrepInfo);
        end
    else
        Urest = speye(dU);
    end
    % extract nontrivial representations
    C = full(U*samples.commutantSample(1)*U');
    [V D] = replab.rep.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    nontrivials = cell(1, n);
    for i = 1:n
        nontrivials{i} = rep.subRepUnitary(V(:, runs{i})' * U, [], replab.IrrepInfo);
    end
    sub = horzcat(trivials, nontrivials);
end
