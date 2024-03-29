function [P, exitFlag] = findProjection_largeScale(rep, I, nSamples, tolerances, Ip, Pp, forceReal)
% Finds a projection map for a subrepresentation defined by an injection map
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   I (double(\*,\*)): Injection map matrix
%   nSamples (integer): Number of samples to use in the large-scale version of the algorithm
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing biorthogonality
%   Pp (double(e,D)): Projection map matrix prescribing biorthogonality
%   forceReal (logical): Whether to force the projection map to be real
%
% Returns:
%   double(\*,\*): Projection map
    D = rep.dimension;
    d = size(I, 2);
    e = size(Ip, 2);
    replab.msg(1, 'Projection map search: dim(parent) = %d, dim(subrep) = %d', D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   ortho    delta');
    replab.msg(2, '-----------------------');
    if forceReal
        P = replab.numerical.randomUnitaryOver(d, 'R') * I'; % as good as a guess as anything else
    else
        P = replab.numerical.randomUnitaryOver(d, rep.field) * I'; % as good as a guess as anything else
    end
    [P, ~] = replab.rep.biorthoStepP(I, P);
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    tolerances.logHeader;
    while exitFlag == 0
        I1 = zeros(D, d);
        P1 = zeros(d, D);
        for j = 1:nSamples
            g = rep.group.sample;
            rhoI = rep.matrixRowAction(g, I);
            Prho = rep.matrixColAction(g, P);
            P1 = P1 + (P * rhoI) * Prho;
        end
        if forceReal
            P1 = real(P1);
        end
        if e > 0
            P1 = P1 - (P1 * Ip) * Pp;
        end
        P1 = P1/nSamples;
        P1 = replab.rep.biorthoStepP(I, P1);
        delta(k) = norm(P1 - P, 'fro')/norm(P, 'fro');
        omega(k) = delta(k);
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        P = P1;
    end
end
