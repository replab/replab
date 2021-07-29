function [v1, v2, exitFlag] = split_unitary_largeScale_old(sub, nSamples, tolerances, Ip)
% Refines a generic unitary subrepresentation
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation to split
%   nSamples (integer): Number of samples to use for each approximate averaging (rec: 5)
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing orthogonality
    rho = sub.parent;
    D = sub.parent.dimension;
    d = sub.dimension;
    e = size(Ip, 2);
    replab.msg(1, 'Unitary split over %s: dim(parent) = %d, dim(subrep) = %d', sub.field, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    assert(rho.isUnitary);
    %if ~gen.mapsAreAdjoint
    %    [I0, ~] = replab.numerical.econqr(I0);
    %end
    k = 1;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    tolerances.logHeader;
    v1 = randn(d, 1);
    v2 = randn(d, 1);
    if sub.overC
        v1 = v1 + 1i * randn(d, 1);
        v2 = v2 + 1i * randn(d, 1);
    end
    v1 = v1/norm(v1);
    v2 = v2/norm(v2);
    v1 = sub.injection * v1;
    v2 = sub.injection * v2;
    while exitFlag == 0
        x = randn(d, 1);
        if sub.overC
            x = x + 1i * randn(d, 1);
        end
        x = x/norm(x);
        x = sub.injection * x;
        x1 = zeros(D, 1);
        for j = 1:nSamples
            g = rho.group.sample;
            rhov1 = rho.matrixRowAction(g, v1);
            rhov2 = rho.matrixRowAction(g, v2);
            x1 = x1 + rhov1 * (rhov1' * x) + rhov2 * (rhov2' * x);
        end
        x1 = x1/norm(x1);
        x1 = x;

        g = cell(1, nSamples);
        for j = 1:nSamples
            g{j} = rho.group.sample;
        end

        w1p = zeros(D, 1);
        for j = 1:nSamples
            rhov1 = rho.matrixRowAction(g{j}, v1);
            rhov2 = rho.matrixRowAction(g{j}, v2);
            w1p = w1p + rhov1 * (rhov1' * x1) - rhov2 * (rhov2' * x1);
        end
        w1p = w1p/nSamples;

        alpha1 = w1p'*x1;
        w1 = w1p - alpha1*x1;

        beta2 = norm(w1);
        x2 = w1/beta2;

        w2p = zeros(D, 1);
        for j = 1:nSamples
            rhov1 = rho.matrixRowAction(g{j}, v1);
            rhov2 = rho.matrixRowAction(g{j}, v2);
            w2p = w2p + rhov1 * (rhov1' * x2) - rhov2 * (rhov2' * x2);
        end
        w2p = w2p/nSamples;

        alpha2 = w2p'*x2;
        w2 = w2p - alpha2*x2 - beta2*x1;

        V = [x1 x2];
        T = [alpha1 beta2
             beta2 alpha2];
        [U, ev] = eig(T);
        newv = V*U;
        v1 = newv(:,1);
        v2 = newv(:,2);
        if e > 0
            v1 = v1 - Ip * (Ip' * v1);
        end
        if e > 0
            v2 = v2 - Ip * (Ip' * v2);
        end
        v1 = v1/norm(v1);
        v2 = v2/norm(v2);
        delta(k) = abs(trace(ev));
        omega(k) = norm(w2);
        exitFlag = tolerances.test(omega, delta, k);
        % Force I1 to be unitary using a single Newton iteration
        k = k + 1;
    end
end
