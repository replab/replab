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
            [I, ~] = qr(I, 0); % find orthogonal basis of the range of I
            P = I';
        end
    else
        I = sub.injection;
        P = sub.projection;
    end
    % Get a sample from the parent representation commutant
    [X, projErr] = replab.irreducible.sampleCommutantRealEV(rep);
    % Sample in the subrepresentation
    A = P * X * I;
    % Find the eigendecomposition, special path for unitary parent rep.
    if unitary
        A = (A + A')/2;
        [V, D] = eig(A, 'vector');
        W = V';
    else
        [V, D, W] = replab.numerical.realeig(A);
        D = full(diag(D));
        W = W';
    end
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
    end
    % Identify clusters of eigenvalues
    start = 1;
    blocks = cell(1, 0);
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
        blocks{1,end+1} = start:(next-1);
        start = next;
    end
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
