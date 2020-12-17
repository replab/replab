function subs = split(rep, sub, safetyFactor)
% Splits a representation into subrepresentations by using a commutant sample
%
% Decomposes a possibly non-irreducible subrepresentation ``rep`` into subrepresentations.
%
% Args:
%   rep (`+replab.Rep`): Representation of which to decompose a subspace
%   sub (`+replab.SubRep`): Subspace to decompose further
%   safetyFactor (double, optional): Factor on the error estimation using eigendecomposition residuals, default: 100
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Subrepresentations in ``sub`` with their ``.parent`` set to ``rep``
    if nargin < 3 || isempty(safetyFactor)
        safetyFactor = 100;
    end
    replab.log(1, '*** Split subspace using commutant sample: dim(parent) = %d, dim(subspace) = %d', rep.dimension, sub.dimension);
    d0 = rep.dimension;
    d1 = sub.dimension;
    I = sub.injection;
    P = sub.projection;
    % Get a sample from the parent representation commutant
    % Sample from the subrepresentation space
    dom = replab.domain.Matrices(sub.field, sub.dimension, sub.dimension);
    X = I * dom.sample * P;
    t = cputime;
    [X, projErr] = rep.commutant.project(X);
    % Sample in the subrepresentation
    A = P * X * I;
    % Find the eigendecomposition
    t = cputime;
    [V, D, W] = eig(A, 'vector');
    blocks = cell(1, 0);
    values = zeros(1, 0);
    i = 1;
    while i <= d1
        values(1, end+1) = real(D(i));
        if i < d1 && real(D(i)) == real(D(i+1))
            blocks{1, end+1} = [i i+1];
            i = i + 2;
        else
            blocks{1, end+1} = i;
            i = i + 1;
        end
    end
    [~, I] = sort(values, 'ascend');
    blocks = blocks(I);
    ind = [blocks{:}];
    V = sub.injection * V(:,ind);
    W = W(:,ind)' * sub.projection;
    D = D(ind);
    % Compute residuals according to http://www.netlib.org/utk/people/JackDongarra/etemplates/node278.html
    R = X * V - V * diag(sparse(D)); % column vectors are residuals for right eigenvectors
    S = W * X - diag(sparse(D)) * W; % row vectors are residuals for left eigenvectors
    evError = zeros(1, d1); % Error estimates
    for i = 1:d1
        evError(i) = max(norm(S(i,:)), norm(R(:,i)))/abs(W(i,:)*V(:,i));
        evError(i) = max(evError(i), eps(D(i)));
    end
    replab.log(2, 'Estimated error on eigenvalues: min %e mean %e max %e', min(evError), mean(evError), max(evError));
    % Identify clusters of eigenvalues
    start = 1;
    blocks = cell(1, 0);
    blockError = zeros(1, 0);
    while start <= d1
        next = start + 1;
        maxEVE = evError(start);
        % Identify close eigenvalues forwards
        while next <= d1 && real(D(next) - D(next-1)) < safetyFactor*(maxEVE + max(maxEVE, evError(next))) + projErr
            maxEVE = max(maxEVE, evError(next));
            next = next + 1;
        end
        % Identify close eigenvalues backwards, in case the max. error estimate of EV so far is bigger
        while start > 1 && real(D(start) - D(start-1)) < safetyFactor*(maxEVE + max(maxEVE, evError(start-1))) + projErr
            lastBlock = blocks{end};
            maxEVE = max(maxEVE, evError(lastBlock));
            start = lastBlock(1);
            blocks = blocks(1:end-1);
        end
        block = start:(next-1);
        blocks{1,end+1} = block;
        blockError(1,end+1) = max(real(D(block))) - min(real(D(block)));
        start = next;
    end
    replab.log(1, 'Identified %d invariant subspaces', length(blocks));
    replab.log(2, 'Eigenvalue spread in clusters: min %e mean %e max %e', min(blockError), mean(blockError), max(blockError));
    % Create subrepresentations
    subs = cell(1, length(blocks));
    for i = 1:length(blocks)
        block = blocks{i};
        injection = V(:, block);
        projection = W(block, :);
        assert(isreal(D(block)) == isreal(injection));
        assert(isreal(D(block)) == isreal(projection));
        if ~isreal(D(block))
            i1 = injection(:, 1:2:end);
            i2 = injection(:, 2:2:end);
            p1 = projection(1:2:end, :);
            p2 = projection(2:2:end, :);
            i3 = i1 + i2;
            p3 = p1 + p2;
            i4 = 1i*(i1 - i2);
            p4 = 1i*(p1 - p2);
            assert(isreal(i3) && isreal(i4) && isreal(p3) && isreal(p4));
            injection = [i3 i4];
            projection = [p3; p4];
            projection = (projection*injection)\projection;
            encodesComplexStructure = true;
        else
            encodesComplexStructure = false;
        end
        subs{i} = rep.subRep(injection, 'projection', projection, 'encodesComplexStructure', encodesComplexStructure);
    end
end
