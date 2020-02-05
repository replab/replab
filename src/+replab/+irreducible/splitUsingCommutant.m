function sub = splitUsingCommutant(rep, samples, sub)
% splitUsingCommutant - Splits a representation into irreducible representations by using a commutant sample
%
% ``sub = replab.irreducible.splitUsingCommutant(rep, samples, sub)``
%
% Decomposes a possibly non-irreducible subrepresentation of ``rep`` into irreducible subrepresentations.
% Does not group the irreducible subrepresentations into isotypic components, but identifies every instance of
% the trivial representation by setting ``label`` to ``'1'`` in the subrepresentation `+replab.IrrepInfo`.
%
% Args:
%   rep (`+replab.Rep`): Parent representation to decompose
%   samples (`+replab.+irreducible.OnDemandSamples`): Samples from the parent representation
%   sub (`+replab.SubRep`): Subrepresentation of ``rep`` to decompose
%
% Returns:
%   row cell array of `+replab.SubRep`: Irreducible subrepresentations inside ``sub`` as subrepresentations of ``rep``
%
% Notes:
%   We first extract the trivial component, as the irreducible construction depends on having
%   that component identified. We then split the orthogonal subspace using a commutant sample.
%   At each step of the process below, we attempt rational basis recovery and construct
%   subrepresentations accordingly so that "nice bases" are presented to user as much as
%   possible.
    d = rep.dimension;
    replab.irreducible.tell('splitUsingCommutant dimension %d', sub.dimension);
    dSub = sub.dimension;
    tol = replab.Parameters.doubleEigTol;
    trivialIrrepInfo = replab.irreducible.TrivialInfo(rep.field);
    nontrivialIrrepInfo = replab.irreducible.Info([], []);
    trivials = {};
    % extract trivial representations
    S = sub.U*samples.trivialSample(1)*sub.U';
    if replab.isNonZeroMatrix(S, replab.Parameters.doubleEigTol) % magic epsilon
        UTrivial = orth(S')';
        nTrivial = size(UTrivial, 1);
        trivial = cell(1, nTrivial);
        if rep.overC
            % try to recover real coefficients for the trivial basis
            rec = replab.irreducible.recoverReal(UTrivial);
            if ~isempty(rec)
                UTrivial = rec;
            end
        end
        % try to recover integer coefficients for the trivial basis
        assert(size(UTrivial, 2) == dSub);
        if isreal(UTrivial)
            [VTrivial isOrtho] = replab.rational.integerRowSpan(UTrivial);
            assert(size(VTrivial, 1) == nTrivial);
        else
            VTrivial = [];
        end
        if ~isempty(VTrivial) && isOrtho
            % success: preserve integer coefficients if they exist
            for i = 1:nTrivial
                trivials{i} = sub.subRepUnitaryByIntegerBasis(VTrivial(i,:), trivialIrrepInfo).collapseParent;
            end
        else
            % failure: just use the unitary basis we already have
            for i = 1:nTrivial
                trivials{i} = sub.subRepUnitary(UTrivial(i,:), [], trivialIrrepInfo).collapseParent;
            end
        end
        % we still need to decompose the orthogonal complement
        URest = null(UTrivial)';
        nNonTrivial = size(URest, 1);
        assert(size(URest, 2) == dSub);
        assert(nTrivial + nNonTrivial == dSub);
    else
        % no trivial component? what remains to decompose is the whole subspace
        URest = speye(dSub);
    end
    % extract nontrivial representations by sampling the commutant
    % however, we use a sample of the parent representation, which can be computed often faster,
    % and shared among computations
    C = full(URest*sub.U*samples.commutantSample(1)*sub.U'*URest');
    C = C+C';
    % the eigenspace decomposition is the basis of the numerical decomposition
    [UMat D] = replab.irreducible.sortedEig(C, 'ascend', false);
    D = diag(D);
    D = D(:)';
    mask = bsxfun(@(x,y) abs(x-y)<tol, D, D');
    % find repeated eigenvalues
    runs = replab.Partition.connectedComponents(mask).blocks;
    n = length(runs);
    nontrivials = cell(1, n);
    for i = 1:n
        UIrrep = UMat(:, runs{i})' * URest;
        if rep.overC
            % try to recover real coefficients for the trivial basis
            rec = replab.irreducible.recoverReal(UIrrep);
            if ~isempty(rec)
                UIrrep = rec;
            end
        end
        % try to recover an integer basis
        if isreal(UIrrep)
            VIrrep = replab.rational.integerRowSpan(UIrrep);
        else
            VIrrep = [];
        end
        if ~isempty(VIrrep)
            nontrivials{i} = sub.subRepUnitaryByIntegerBasis(VIrrep, nontrivialIrrepInfo).collapseParent;
        else
            nontrivials{i} = sub.subRepUnitary(UIrrep, [], nontrivialIrrepInfo).collapseParent;
        end
    end
    sub = horzcat(trivials, nontrivials);
end
