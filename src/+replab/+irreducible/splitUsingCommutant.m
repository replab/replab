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
    d = rep.dimension;
    knownUnitary = isequal(rep.isUnitary, true);
    if ~knownUnitary
        [A Ainv] = rep.unitaryChangeOfBasis;
        % A: unitary space <- nonunitary space
        % Ainv: nonunitary space <- unitary space
    else
        rep1 = rep;
        A = speye(d);
        Ainv = speye(d);
    end
    tol = replab.Parameters.doubleEigTol;
    % extract nontrivial representations by sampling the commutant
    C1 = A * rep.commutant.sampleInContext(context, 1) * Ainv;
    % the eigenspace decomposition is the basis of the numerical decomposition
    % V'*C*V = D
    [U1 D] = replab.irreducible.sortedEig((C1 + C1')/2, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    % find repeated eigenvalues
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    if n == 1
        rep.isIrreducible = true;
    end
    sub = cell(1, n);
    for i = 1:n
        basis = Ainv * U1(:, runs{i});
        embedding = U1(:, runs{i})' * A;
        sub{1,i} = rep.subRep(basis, embedding);
        sub{1,i}.isIrreducible = true;
    end
end
