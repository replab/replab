function gen = harmonize_unitary_largeScale(gen0, genM, nSamples, tolerances, Ip)
% Refines a generic non-unitary subrepresentation
%
% Example:
%   >>> G = replab.S(3);
%   >>> genM = replab.rep.GenSubRep.fromSubRep(G.standardRep.unitarize.collapse);
%   >>> gen0 = replab.rep.GenSubRep.fromSubRep(G.standardRep.withNoise(0.1, 0.1));
%   >>> gen = replab.rep.harmonize_unitary_largeScale(gen0, genM, 5, replab.rep.Tolerances, []);
%   >>> g = [3 2 1];
%   >>> norm(genM.image(g) - gen.image(g)) < 1e-10
%       1
%
% Args:
%   gen0 (`+replab.+rep.GenSubRep`): Generic subrepresentation to harmonize
%   genM (`+replab.+rep.GenSubRep`): Model generic subrepresentation
%   nSamples (integer): Number of samples per averaging iteration
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing orthogonality
%
% Returns:
%   `+replab.+rep.GenSubRep`: Refined generic subrepresentation
    D = gen0.parent.dimension;
    d = gen0.dimension;
    rho = gen0.parent;
    mu = genM.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Unitary harmonization over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    tolerances.logHeader;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    Im = genM.injection;
    I = gen0.injection * replab.numerical.randomUnitaryOver(d, gen0.divisionRing);
    while exitFlag == 0
        I1 = zeros(D, d);
        for j = 1:nSamples
            g = rho.group.sample;
            switch type
              case {'R/R', 'C/C'}
                muIm = mu.matrixRowAction(g, Im);
                rhoI = rho.matrixRowAction(g, I);
              case 'C/R'
                muIm = mu.matrixRowAction(g, real(Im)) + ...
                       mu.matrixRowAction(g, imag(Im)) * 1i;
                rhoI = rho.matrixRowAction(g, real(I)) + ...
                       rho.matrixRowAction(g, imag(I)) * 1i;
              case 'H/R'
                muIm = replab.H(mu.matrixRowAction(g, part1(Im)), ...
                                mu.matrixRowAction(g, parti(Im)), ...
                                mu.matrixRowAction(g, partj(Im)), ...
                                mu.matrixRowAction(g, partk(Im)));
                rhoI = replab.H(rho.matrixRowAction(g, part1(I)), ...
                                rho.matrixRowAction(g, parti(I)), ...
                                rho.matrixRowAction(g, partj(I)), ...
                                rho.matrixRowAction(g, partk(I)));
            end
            I1 = I1 + rhoI * (muIm' * Im);
        end
        if ~isempty(Ip)
            I1 = I1 - Ip * (Ip' * I1);
        end
        f = trace(I1*I1')/d;
        I1 = I1/sqrt(f);
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        omega(k) = norm(I1'*I1 - eye(d), 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
    end
    % One step of Newton iteration to force unitarity (helps!)
    N = I'*I;
    P = I*N/2;
    I = 2*I + P*N - 3*P;
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, true, I, I');
end
