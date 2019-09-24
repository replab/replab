function sub = splitUsingCommutant(rep, samples, U)
% Splits a representation into irreducible representations by using a commutant sample
    d = rep.dimension;
    if nargin < 3
        U = speye(d);
    end
    replab.irreducible.tell('splitUsingCommutant dimension %d', size(U, 1));
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
    if replab.isNonZeroMatrix(S, replab.Settings.doubleEigTol) % magic epsilon
        Utrivial = orth(S')';
        % try to recover real indices
        rec = replab.irreducible.recoverReal(Utrivial);
        if ~isempty(rec)
            Utrivial = rec;
        end
        Urest = null(Utrivial)';
        assert(size(Utrivial, 2) == dU);
        assert(size(Urest, 2) == dU);
        assert(size(Urest, 1) + size(Utrivial, 1) == dU);
        for i = 1:size(Utrivial, 1)
            trivials{1,i} = rep.subRepUnitary(Utrivial(i,:) * U, [], trivialIrrepInfo);
        end
    else
        Urest = speye(dU);
    end
    % extract nontrivial representations
    C = full(Urest*U*samples.commutantSample(1)*U'*Urest');
    C = C+C';
    [V D] = replab.irreducible.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    nontrivials = cell(1, n);
    for i = 1:n
        nontrivials{i} = rep.subRepUnitary(V(:, runs{i})' * Urest * U, [], replab.IrrepInfo);
    end
    sub = horzcat(trivials, nontrivials);
end
