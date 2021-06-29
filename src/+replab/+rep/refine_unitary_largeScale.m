function [gen, exitFlag] = refine_unitary_largeScale(gen0, nSamples, tolerances, Ip)
% Refines a generic unitary subrepresentation
%
% Example:
%   >>> G = replab.S(3);
%   >>> rep = G.standardRep.unitarize.collapse.withNoise(0.1);
%   >>> rep.projectorErrorBound > 1e-3
%       1
%   >>> gen1 = replab.rep.refine_unitary_largeScale(replab.rep.GenSubRep.fromSubRep(rep), 5, replab.rep.Tolerances, []);
%   >>> gen1.toSubRep.projectorErrorBound < 1e-10
%       1
%
% Args:
%   gen0 (`+replab.+rep.GenSubRep`): Generic subrepresentation to refine
%   nSamples (integer): Number of samples per averaging iteration
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing orthogonality
%
% Returns:
%   `+replab.+rep.GenSubRep`: Refined generic subrepresentation
    D = gen0.parent.dimension;
    d = gen0.dimension;
    e = size(Ip, 2);
    rho = gen0.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Unitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    I0 = gen0.injection;
    I = I0;
    assert(rho.isUnitary);
    %if ~gen.mapsAreAdjoint
    %    [I0, ~] = replab.numerical.econqr(I0);
    %end
    k = 1;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    tolerances.logHeader;
    while exitFlag == 0
        I1 = zeros(size(I));
        for j = 1:nSamples
            g = rho.group.sample;
            switch type
              case {'R/R', 'C/C'}
                rhoI = rho.matrixRowAction(g, I);
              case 'C/R'
                rhoI = rho.matrixRowAction(g, real(I)) + ...
                       rho.matrixRowAction(g, imag(I)) * 1i;
              case 'H/R'
                rhoI = replab.H(rho.matrixRowAction(g, part1(I)), ...
                                rho.matrixRowAction(g, parti(I)), ...
                                rho.matrixRowAction(g, partj(I)), ...
                                rho.matrixRowAction(g, partk(I)));
            end
            I1 = I1 + rhoI * (rhoI' * I);
        end
        I1 = I1/nSamples;
        if e > 0
            I1 = I1 - Ip * (Ip' * I1);
        end
        f = trace(I1*I1')/d;
        I1 = I1/sqrt(f);
        delta(k) = norm(I1 - I, 'fro')/norm(I1, 'fro');
        omega(k) = norm(I1'*I1 - eye(d), 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        % Force I1 to be orthogonal using a single Newton iteration
        N = I1'*I1;
        P = I1*N/2;
        I1 = 2*I1 + P*N - 3*P;
        I = I1;
        k = k + 1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, true, I, I');
end
