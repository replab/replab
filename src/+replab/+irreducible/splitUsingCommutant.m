function sub = splitUsingCommutant(rep, context)
% Splits a representation into irreducible representations by using a commutant sample
%
% Decomposes a possibly non-irreducible representation ``rep`` into irreducible subrepresentations.
%
% The trivial subrepresentation must already have been taken out.
    if ~isequal(rep.trivialDimension, 0)
        sub = replab.DispatchNext('Need to have no trivial subrepresentation');
        return
    end
    if ~isequal(rep.isUnitary, true)
        sub1 = replab.irreducible.splitUsingCommutant(rep.unitarize);
        if isa(sub1, 'replab.DispatchNext')
            sub = sub1;
            return
        end
        sub = cellfun(@(s) replab.rep.collapse(s), sub1, 'uniform', 0);
        return
    end
    d = rep.dimension;
    tol = replab.Parameters.doubleEigTol;
    % extract nontrivial representations by sampling the commutant
    C = rep.commutant.sampleInContext(context, 1);
    % the eigenspace decomposition is the basis of the numerical decomposition
    % V'*C*V = D
    [U1 D] = replab.numerical.sortedEig((C + C')/2, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    % find repeated eigenvalues
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    if n == 1
        rep.isIrreducible = true;
        sub = {replab.SubRep.fullSubRep(rep)};
        return
    end
    sub = cell(1, n);
    for i = 1:n
        basis = U1(:, runs{i});
        sub{1,i} = rep.subRep(basis, basis');
        sub{1,i}.isIrreducible = true;
    end
end
