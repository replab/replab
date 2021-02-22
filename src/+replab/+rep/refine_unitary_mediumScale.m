function gen = refine_unitary_mediumScale(gen0, nInnerIterations, tolerances, Ip)
% Refines a generic unitary representation
%
% Example:
%   >>> G = replab.S(3);
%   >>> rep = G.standardRep.unitarize.collapse.withNoise(0.1);
%   >>> rep.projectorErrorBound > 1e-3
%       1
%   >>> gen1 = replab.rep.refine_unitary_mediumScale(replab.rep.GenSubRep.fromSubRep(rep), 5, replab.rep.Tolerances, []);
%   >>> gen1.toSubRep.projectorErrorBound < 1e-10
%       1
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to refine
%   nInnerIterations (integer): Number of inner iterations
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing orthogonality
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    D = gen0.parent.dimension;
    d = gen0.dimension;
    rho = gen0.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Unitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Full commutant projection algorithm');
    replab.msg(1, '');
    tolerances.logHeader;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    I0 = gen0.injection;
    if ~gen0.mapsAreAdjoint
        [I0, ~] = replab.numerical.qr(I0);
    end
    I = I0;
    iter = 1;
    while exitFlag == 0
        Ftilde = I*I';
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
        for j = 1:nInnerIterations
            I1 = Foverline * I;
            if ~isempty(Ip)
                I1 = I1 - Ip * (Ip' * I1);
            end
            [I1, ~] = replab.numerical.qr(I1);
        end
        omega(k) = norm(Ftilde - Foverline, 'fro');
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, true, I, I');
end
