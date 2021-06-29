function [I, ev, exitFlag] = split_unitary_largeScale(rho, d, nSamples, tolerances)
% Refines a generic unitary subrepresentation
%
% Args:
%   rho (`+replab.Rep`): Representation to split
%   d (integer): Block size to use, must be ``< rho.dimension``
%   nSamples (integer): Number of samples to use for each approximate averaging (rec: 5)
%   tolerances (`.Tolerances`): Termination criteria
    D = rho.dimension;
    replab.msg(1, 'Unitary split over %s: dim = %d', rho.field, D);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    assert(rho.isUnitary);
    k = 1;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    tolerances.logHeader;
    e = 1; % number of extra dimensions
    I = randn(D, d);
    X = randn(d+e, d+e);
    if rho.overC
        I = I + 1i * randn(D, d);
        X = X + 1i * randn(d+e, d+e);
    end
    X = X + X';
    X = X/trace(X);
    [I, ~] = replab.numerical.econqr(I);
    while exitFlag == 0
        % sample extra basis elements
        assert(e == 1);
        x = randn(D, 1);
        if rho.overC
            x = x + 1i * randn(D, 1);
        end
        x = x - I * (I' * x);
        x = x/norm(x);
        I = [I x];
        X(:,d+1:end) = 0;
        X(d+1:end,:) = 0;
        X1 = X;
        for l = 1:5
            X2 = zeros(d+1, d+1);
            for j = 1:nSamples
                g = rho.group.sample;
                img = I'*rho.matrixRowAction(g, I);
                X2 = X2 + img*X1*img';
            end
            X2 = X2/nSamples;
            X1 = (X2+X2')/2;
        end
        [U, ev] = eig(X1);
        ev = diag(ev);
        [~, ind] = sort(ev, 'desc');
        ev = ev(ind);
        ev(:).'
        U = U(:,ind);
        X1 = U'*X1*U;
        I = I * U(:,1:d);
        X1 = X1/trace(X1);
        delta(k) = norm(X1 - X, 'fro');
        omega(k) = min(ev);
        exitFlag = tolerances.test(omega, delta, k);
        X = X1;
        k = k + 1;
    end
end
