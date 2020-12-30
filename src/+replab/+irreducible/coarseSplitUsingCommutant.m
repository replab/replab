function subs = coarseSplitUsingCommutant(rep, sub, safetyFactor)
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
    replab.msg(1, '*** Split subspace using commutant sample: dim(parent) = %d, dim(subspace) = %d', rep.dimension, sub.dimension);
    d0 = rep.dimension;
    d1 = sub.dimension;
    % If the parent representation is unitary, then use a unitary variant, and make the injection map
    % of the subrepresentation an isometry
    unitary = rep.knownUnitary;
    if unitary
        I = sub.injection;
        if sub.knownUnitary
            P = I';
        else
            replab.msg(2, 'Subrepresentation is not unitary, but parent representation is: making subrepresentation unitary');
            [I, ~] = qr(I, 0); % find orthogonal basis of the range of I
            P = I';
        end
    else
        I = sub.injection;
        P = sub.projection;
    end
    % Get a sample from the parent representation commutant

    % sample from the subrepresentation space
    dom = replab.domain.Matrices(sub.field, sub.dimension, sub.dimension);
    X = I * dom.sample * P;
    if ~rep.knownUnitary
        % make eigenvalues real in the basis where such matrices are self-adjoint
        U = rep.unitarize;
        X = U.A * X * U.Ainv;
        X = (X + X')/2;
        X = U.Ainv * X * U.A;
    else
        X = (X + X')/2;
    end
    t = cputime;
    [X, projErr] = rep.commutant.project(X);
    replab.msg(2, 'Time (commutant projection): %2.2f s', cputime - t);
    replab.msg(2, 'Error upper bound on the projection (Frob. norm): %e', projErr);
    % Sample in the subrepresentation
    A = P * X * I;
    % Find the eigendecomposition, special path for unitary parent rep.
    t = cputime;
    if unitary
        A = (A + A')/2;
        [V, D] = eig(A, 'vector');
        W = V';
    else
        [V, D, W] = replab.numerical.realeig(A);
        D = full(diag(D));
        W = W';
    end
    replab.msg(2, 'Time (eigendecomposition): %2.2f s', cputime - t);
    % Sort eigenvalues
    D = D(:)';
    [~, I] = sort(D, 'ascend');
    V = sub.injection * V(:,I);
    W = W(I,:) * sub.projection;
    D = D(I);
    % Compute residuals according to http://www.netlib.org/utk/people/JackDongarra/etemplates/node278.html
    R = X * V - V * diag(sparse(D)); % column vectors are residuals for right eigenvectors
    S = W * X - diag(sparse(D)) * W; % row vectors are residuals for left eigenvectors
    evError = zeros(1, d1); % Error estimates
    for i = 1:d1
        evError(i) = max(norm(S(i,:)), norm(R(:,i)))/abs(W(i,:)*V(:,i));
        evError(i) = max(evError(i), eps(D(i)));
    end
    replab.msg(2, 'Estimated error on eigenvalues: min %e mean %e max %e', min(evError), mean(evError), max(evError));
    % Identify clusters of eigenvalues
    start = 1;
    blocks = cell(1, 0);
    blockError = zeros(1, 0);
    while start <= d1
        next = start + 1;
        maxEVE = evError(start);
        % Identify close eigenvalues forwards
        while next <= d1 && (D(next) - D(next-1)) < safetyFactor*(maxEVE + max(maxEVE, evError(next))) + projErr
            maxEVE = max(maxEVE, evError(next));
            next = next + 1;
        end
        % Identify close eigenvalues backwards, in case the max. error estimate of EV so far is bigger
        while start > 1 && D(start) - D(start-1) < safetyFactor*(maxEVE + max(maxEVE, evError(start-1))) + projErr
            lastBlock = blocks{end};
            maxEVE = max(maxEVE, evError(lastBlock));
            start = lastBlock(1);
            blocks = blocks(1:end-1);
        end
        block = start:(next-1);
        blocks{1,end+1} = block;
        blockError(1,end+1) = max(D(block)) - min(D(block));
        start = next;
    end
    replab.msg(1, 'Identified %d invariant subspaces', length(blocks));
    replab.msg(2, 'Eigenvalue spread in clusters: min %e mean %e max %e', min(blockError), mean(blockError), max(blockError));
    % Create subrepresentations
    subs = cell(1, length(blocks));
    for i = 1:length(blocks)
        block = blocks{i};
        injection = V(:, block);
        if unitary
            subs{i} = rep.subRep(injection, 'projection', injection', 'isUnitary', true);
        else
            projection = W(block, :);
            projection = (projection*injection)\projection;
            subs{i} = rep.subRep(injection, 'projection', projection);
        end
    end
end
