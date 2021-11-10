function gen = refine_nonUnitary_mediumScale(gen0, nInnerIterations, tolerances, Ip, Pp)
% Refines a generic non-unitary representation
%
% Example:
% %   >>> G = replab.S(3);
% %   >>> rep = G.standardRep.withNoise(0.1);
% %   >>> rep.projectorErrorBound > 1e-3
% %       1
% %   >>> gen1 = replab.rep.refine_nonUnitary_mediumScale(replab.rep.GenSubRep.fromSubRep(rep), 5, replab.rep.Tolerances, []);
% %   >>> gen1.toSubRep.projectorErrorBound < 1e-10
% %       1
%
% Args:
%   gen0 (`+replab.+rep.GenSubRep`): Generic subrepresentation to refine
%   nInnerIterations (integer): Number of inner iterations
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
    replab.msg(1, 'Using the full commutant projection algorithm');
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
        Ftilde = I*P;
        switch type
          case {'R/R', 'C/C'}
            [Foverline, err] = rho.commutant.project(Ftilde);
          case 'C/R'
            [F1, e1] = rho.commutant.project(real(Ftilde));
            [Fi, ei] = rho.commutant.project(imag(Ftilde));
            Foverline = F1 + Fi*1i;
            err = sqrt(e1^2 + ei^2);
          case 'H/R'
            [F1, e1] = rho.commutant.project(Ftilde.part1);
            [Fi, ei] = rho.commutant.project(Ftilde.parti);
            [Fj, ej] = rho.commutant.project(Ftilde.partj);
            [Fk, ek] = rho.commutant.project(Ftilde.partk);
            Foverline = replab.H(F1, Fi, Fj, Fk);
            err = sqrt(e1^2 + ei^2 + ej^2 + ek^2);
        end
        I1 = I;
        P1 = P;
        for j = 1:nInnerIterations
            I1 = Foverline * I1;
            P1 = P1 * Foverline;
            if ~isempty(Ip) && ~isempty(Pp)
                I1 = I1 - Ip * (Pp * I1);
                P1 = P1 - (P1 * Ip) * Pp;
            end
            I1 = I1 / (P0 * I1);
            P1 = (P1 * I1) \ P1;
        end
        omega(k) = norm(Ftilde - Foverline, 'fro');
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
        P = P1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, false, I, P);
end
