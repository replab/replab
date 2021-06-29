function gen = harmonize_nonUnitary_largeScale(gen0, genM, nSamples, tolerances, Ip, Pp)
% Refines a generic non-unitary subrepresentation
%
% Example:
%   >>> G = replab.S(3);
%   >>> genM = replab.rep.GenSubRep.fromSubRep(G.standardRep);
%   >>> gen0 = replab.rep.GenSubRep.fromSubRep(G.standardRep.unitarize.collapse.withNoise(0.1));
%   >>> gen = replab.rep.harmonize_nonUnitary_largeScale(gen0, genM, 5, replab.rep.Tolerances, [], []);
%   >>> g = [3 2 1];
%   >>> norm(genM.image(g) - gen.image(g)) < 1e-10
%       1
%
% Args:
%   gen (`+replab.+rep.GenSubRep`): Generic subrepresentation to harmonize
%   genM (`+replab.+rep.GenSubRep`): Model generic subrepresentation
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
    mu = genM.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Nonunitary harmonization over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    tolerances.logHeader;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    Im = genM.injection;
    Pm = genM.projection;
    I = gen0.injection;
    P = gen0.projection;
    while exitFlag == 0
        I1 = zeros(D, d);
        for j = 1:nSamples
            g = rho.group.sample;
            switch type
              case {'R/R', 'C/C'}
                Pmmu = mu.matrixColAction(g, Pm);
                rhoI = rho.matrixRowAction(g, I);
              case 'C/R'
                Pmmu = mu.matrixColAction(g, real(Pm)) + ...
                       mu.matrixColAction(g, imag(Pm)) * 1i;
                rhoI = rho.matrixRowAction(g, real(I)) + ...
                       rho.matrixRowAction(g, imag(I)) * 1i;
              case 'H/R'
                Pmmu = replab.H(mu.matrixColAction(g, part1(Pm)), ...
                                mu.matrixColAction(g, parti(Pm)), ...
                                mu.matrixColAction(g, partj(Pm)), ...
                                mu.matrixColAction(g, partk(Pm)));
                rhoI = replab.H(rho.matrixRowAction(g, part1(I)), ...
                                rho.matrixRowAction(g, parti(I)), ...
                                rho.matrixRowAction(g, partj(I)), ...
                                rho.matrixRowAction(g, partk(I)));
            end
            I1 = I1 + rhoI * (Pmmu * Im);
        end
        P1 = zeros(d, D);
        for j = 1:nSamples
            g = rho.group.sample;
            switch type
              case {'R/R', 'C/C'}
                muIm = mu.matrixRowAction(g, Im);
                Prho = rho.matrixColAction(g, P);
              case 'C/R'
                muIm = mu.matrixRowAction(g, real(Im)) + ...
                       mu.matrixRowAction(g, imag(Im)) * 1i;
                Prho = rho.matrixColAction(g, real(P)) + ...
                       rho.matrixColAction(g, imag(P)) * 1i;
              case 'H/R'
                muIm = replab.H(mu.matrixRowAction(g, part1(Im)), ...
                                mu.matrixRowAction(g, parti(Im)), ...
                                mu.matrixRowAction(g, partj(Im)), ...
                                mu.matrixRowAction(g, partk(Im)));
                Prho = replab.H(rho.matrixColAction(g, part1(P)), ...
                                rho.matrixColAction(g, parti(P)), ...
                                rho.matrixColAction(g, partj(P)), ...
                                rho.matrixColAction(g, partk(P)));
            end
            P1 = P1 + (Pm * muIm) * Prho;
        end
        if ~isempty(Ip) && ~isempty(Pp)
            I1 = I1 - Ip * (Pp * I);
            P1 = P1 - (P * Ip) * Pp;
        end
        f = trace(P1*I1)/d;
        sf = sqrt(abs(f));
        I1 = I1/sf;
        P1 = P1/sf;
        P1 = P1/sign(f);
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        omega(k) = norm(P1*I1 - eye(d), 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
        P = P1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, genM.isUnitary, I, P);
end
