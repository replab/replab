function irreps = splitUsingCommutantOld(rep, sub, samples, ind, failureProb, errorMargin)
% Splits a representation into irreducible representations by using a commutant sample
%
% Decomposes a possibly non-irreducible subrepresentation ``rep`` into irreducible subrepresentations.
%
% Args:
%   rep (`+replab.Rep`): Representation of which to decompose a subspace
%   sub (`+replab.SubRep`): Subspace to decompose further
%   samples (`.RealEVCommutantSamples`): Samples
%   ind (integer): Sample index
%   failureProb (double): Upper bound on the failure probability
%   errorMargin (double): Margin on the error estimation using
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Irreducible subrepresentations in ``sub`` with their ``.parent`` set to ``rep``
    if nargin < 6 || isempty(errorMargin)
        errorMargin = 10;
    end
    if nargin < 5 || isempty(failureProb)
        failureProb = 1e-9;
    end
    if nargin < 4 || isempty(ind)
        ind = 1;
    end
    if nargin < 3 || isempty(samples)
        samples = replab.irreducible.RealEVCommutantSamples(rep);
    end
    d0 = rep.dimension;
    d1 = sub.dimension;
    [S, projErr] = samples.sample(ind);
    A = sub.projection * S * sub.injection;
    if rep.knownUnitary && sub.knownUnitary
        A = (A + A')/2;
        [V, D] = eig(A, 'vector');
        W = V';
    else
        [V, D, W] = replab.numerical.realeig(A);
        D = full(diag(D));
        W = W';
    end
    D = D(:)';
    [~, I] = sort(abs(D), 'descend');
    V = V(:,I);
    W = W(I,:);
    D = D(I);
    R = A * V - V * diag(sparse(D)); % column vectors are residuals for right eigenvectors
    S = W * A - diag(sparse(D)) * W; % row vectors are residuals for left eigenvectors
    evError = zeros(1, d0);
    for i = 1:d0
        evError(i) = max(norm(S(i,:)), norm(R(:,i)))/abs(W(i,:)*V(:,i));
    end
    % Bound on eigenvalue separation from
    % https://math.stackexchange.com/questions/3207232/minimum-distance-between-independent-points-drawn-from-the-normal-distribution
    %
    % rescale by sqrt(d1) to try to cater to the worst case, TODO: revisit this bound
    % remove the null part
    assert(all(abs(D(d1+1:end)) < errorMargin*evError(d1+1:end) + projErr), 'Error estimation failed');
    assert(all(abs(D(1:d1)) > max([evError(d1+1:end) 0]) + projErr), 'Errors too big; impossible to work in the required subspace');
    D = D(1:d1);
    V = V(:,1:d1);
    W = W(1:d1,:);
    evError = evError(1:d1);
    R = R(:,1:d1);
    S = S(1:d1,:);
    [~, I] = sort(D, 'ascend');
    V = V(:,I);
    W = W(I,:);
    R = R(:,I);
    S = S(I,:);
    D = D(I);
    evError = evError(I);
    mu = 2*erfinv(failureProb/nchoosek(d1, 2))/sqrt(d1);
    % identify clusters of eigenvalues
    start = 1;
    blocks = cell(1, 0);
    while start <= d1
        next = start + 1;
        maxEVE = evError(start);
        while next <= d1 && (D(next) - D(next-1)) < errorMargin*(maxEVE + max(maxEVE, evError(next))) + projErr
            maxEVE = max(maxEVE, evError(next));
            next = next + 1;
        end
        while start > 1 && D(start) - D(start-1) < errorMargin*(maxEVE + max(maxEVE, evError(start-1))) + projErr
            lastBlock = blocks{end};
            maxEVE = max(maxEVE, evError(lastBlock));
            start = lastBlock(1);
            blocks = blocks(1:end-1);
        end
        blocks{1,end+1} = start:(next-1);
        start = next;
    end
    % work eigenspace by eigenspace
    irreps = cell(1, 0);
    for i = 1:length(blocks)
        block = blocks{i};
        injection = V(:,block);
        if rep.knownUnitary && sub.knownUnitary
            sub1 = rep.subRep(injection, 'projection', injection', 'isUnitary', true);
        else
            projection = W(block,:);
            projection = (projection*injection)\projection;
            assert(~any(any(isnan(projection))));
            sub1 = rep.subRep(injection, 'projection', projection);
        end
        if sub.inCache('trivialDimension') && sub.trivialDimension == 0
            sub1.cache('trivialDimension', 0, '==');
        end
        if D(block(end)) - D(block(1)) > mu
            mu1 = 2*erfinv(failureProb/nchoosek(sub1.dimension, 2))/sqrt(sub1.dimension);
            D1 = real(eig(sub1.projection * samples.sample(ind+1) * sub1.injection));
            fprintf('Testing, EV delta %e, max error %e\n', D(block(end))-D(block(1)), mu)
            fprintf('         EV delta %e, max error %e\n', max(D1) - min(D1), mu1)
            if max(D1) - min(D1) > mu1
                irreps1 = replab.irreducible.splitUsingCommutant(rep, sub1, samples, ind + 1, failureProb, errorMargin);
                irreps = horzcat(irreps, irreps1);
            else
                irreps{1,end+1} = sub1;
            end
        else
            fprintf('Good, EV delta %e, max error %e\n', D(block(end))-D(block(1)), mu)
            sub1.cache('isUnitary', true, '==');
            irreps{1,end+1} = sub1;
        end
    end
end
