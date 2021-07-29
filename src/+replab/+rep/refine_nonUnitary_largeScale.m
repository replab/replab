function gen = refine_nonUnitary_largeScale(gen0, nSamples, tolerances, Ip, Pp)
% Refines a generic non-unitary subrepresentation
%
% Example:
%   >>> G = replab.S(3);
%   >>> rep = G.standardRep.withNoise(0.1);
%   >>> rep.projectorErrorBound > 1e-3
%       1
%   >>> gen1 = replab.rep.refine_nonUnitary_largeScale(replab.rep.GenSubRep.fromSubRep(rep), 5, replab.rep.Tolerances, []);
%   >>> gen1.toSubRep.projectorErrorBound < 1e-10
%       1
%
% Args:
%   gen0 (`+replab.+rep.GenSubRep`): Generic subrepresentation to refine
%   nSamples (integer): Number of samples per averaging iteration
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing biorthogonality
%   Pp (double(e,D)): Projection map matrix prescribing biorthogonality
%
% Returns:
%   `+replab.+rep.GenSubRep`: Refined generic subrepresentation
    D = gen0.parent.dimension;
    d = gen0.dimension;
    e = size(Ip, 2);
    rho = gen0.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Nonunitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    tolerances.logHeader;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    I0 = gen0.injection;
    P0 = gen0.projection;
    I = I0;
    P = P0;
    while exitFlag == 0
        I1 = zeros(D, d);
        P1 = zeros(d, D);
        for j = 1:nSamples
            g = rho.group.sample;
            switch type
              case {'R/R', 'C/C'}
                Prho = rho.matrixColAction(g, P);
                rhoI = rho.matrixRowAction(g, I);
              case 'C/R'
                Prho = rho.matrixColAction(g, real(P)) + ...
                       rho.matrixColAction(g, imag(P)) * 1i;
                rhoI = rho.matrixRowAction(g, real(I)) + ...
                       rho.matrixRowAction(g, imag(I)) * 1i;
              case 'H/R'
                Prho = rho.matrixColAction(g, part1(P)) + ...
                       rho.matrixColAction(g, parti(P)) * replab.H.i + ...
                       rho.matrixColAction(g, partj(P)) * replab.H.j + ...
                       rho.matrixColAction(g, partk(P)) * replab.H.k;
                rhoI = rho.matrixRowAction(g, part1(I)) + ...
                       rho.matrixRowAction(g, parti(I)) * replab.H.i + ...
                       rho.matrixRowAction(g, partj(I)) * replab.H.j + ...
                       rho.matrixRowAction(g, partk(I)) * replab.H.k;
            end
            I1 = I1 + rhoI * (Prho * I);
            P1 = P1 + (P * rhoI) * Prho;
        end
        if e > 0
            I1 = I1 - Ip * (Pp * I1);
            P1 = P1 - (P1 * Ip) * Pp;
        end
        f = trace(P1*I1)/d;
        I1 = I1/sqrt(abs(f));
        P1 = P1/(f/sqrt(abs(f)));
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        omega(k) = norm(P1*I1 - eye(d), 'fro');
        if omega(k) < 0.9
            I1 = 2*I1 - I1 * P0 * I1;
            P1 = 2*P1 - P1 * I1 * P1;
        else
            I1 = I1 / (P0 * I1);
            P1 = (P1 * I1) \ P1;
        end
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
        P = P1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, false, I, P);
end
