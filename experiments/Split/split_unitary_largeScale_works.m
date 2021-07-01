function [P, exitFlag] = split_unitary_largeScale_works(rho, nSamples, tolerances)
% Refines a generic unitary subrepresentation
%
% Args:
%   rho (`+replab.Rep`): Representation to split
%   nSamples (integer): Number of samples to use for each approximate averaging (rec: 5)
%   tolerances (`.Tolerances`): Termination criteria
    d = rho.dimension;
    replab.msg(1, 'Unitary split over %s: dim = %d', rho.field, d);
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

    v = randn(d, 1);
    if rho.overC
        v = v + 1i * randn(d, 1);
    end
    v = v/norm(v);
    P = v*v';

    while exitFlag == 0
        P1 = zeros(d, d);
        for j = 1:nSamples
            rhog = rho.sample;
            P1 = P1 + rhog * P * rhog';
        end
        P1 = P1/nSamples;
        dP = min(d, trace(P1));
        P1 = P1*P1;
        P1 = P1/trace(P1)*dP;

        v1 = randn(d, 1);
        if rho.overC
            v1 = v1 + 1i * randn(d, 1);
        end
        v1 = v1/norm(v1);

        % Abbreviated initial iteration step
        w1p = P1 * v1;
        alpha1 = w1p'*v1;
        w1 = w1p - alpha1*v1;

        beta2 = norm(w1);
        v2 = w1/beta2;
        w2p = P1 * v2;
        alpha2 = w2p'*v2;
        w2 = w2p - alpha2*v2 - beta2*v1;

        %beta3 = norm(w2);
        %v3 = w2/beta3;
        %w3p = (P1-Q1) * v3;
        %alpha3 = w3p'*v3;
        %w3 = w3p - alpha3*v3 - beta3*v2;

        V = [v1 v2]; % v3];
                     %T = [alpha1 beta2 0
                     %             beta2 alpha2 beta3
                     %             0 beta3 alpha3];
        T = [alpha1 beta2
             beta2 alpha2];
        [U, ev] = eig(T);
        ev = diag(ev);
        if ev(1) < ev(2)
            U = U(:,[2 1]);
            ev = ev([2 1]);
        end
        new = V*U;
        p = new(:,1);
        q = new(:,2);
        P1 = P1 - ev(2)*q*q';
        P1 = (P1+P1')/2;
        P1 = P1 / (p'*P1*p);
        delta(k) = norm(P1 - P, 'fro');
        omega(k) = abs(ev(1) - 1) + abs(ev(2));
        fprintf('%d ', dP);
        exitFlag = tolerances.test(omega, delta, k);
        % Force I1 to be unitary using a single Newton iteration
        k = k + 1;
        P = P1;
    end
end
