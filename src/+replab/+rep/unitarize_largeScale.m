function gen = unitarize_largeScale(gen0, nSamples, tolerances, Ip, Pp)
% Refines a generic non-unitary subrepresentation
%
% Example:
%   >>> G = replab.S(7);
%   >>> gen0 = replab.rep.GenSubRep.fromSubRep(G.standardRep);
%   >>> gen = replab.rep.unitarize_largeScale(gen0, 5, replab.rep.Tolerances, [], []);
%   >>> g = G.sample;
%   >>> norm(gen.image(g) - gen.inverseImage(g)') < 1e-10
%       1
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to unitarize
%   nSamples (integer): Number of samples per averaging iteration
%   tolerances (`.Tolerances`): Termination criteria
%   Ip (double(D,e)): Injection map matrix prescribing biorthogonality
%   Pp (double(e,D)): Projection map matrix prescribing biorthogonality
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    D = gen0.parent.dimension;
    d = gen0.dimension;
    rho = gen0.parent;
    type = [gen0.divisionRing '/' rho.field];
    replab.msg(1, 'Unitarization over %s: dim(parent) = %d, dim(subrep) = %d', gen0.divisionRing, D, d);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    tolerances.logHeader;
    delta = zeros(1, tolerances.maxIterations);
    omega = zeros(1, tolerances.maxIterations);
    exitFlag = 0;
    k = 1;
    I = gen0.injection;
    P = gen0.projection;
    while exitFlag == 0
        X = zeros(d, d);
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
            img = P*rhoI;
            X = X + img*img';
        end
        X = X/trace(X)*d;
        % approximation of the matrix square root (see Taylor expansion of sqrt(1+x))
        % S = (eye(dsub) + X)/2;
        % and the square root inverse (see Taylor expansion of 1/sqrt(1+x))
        % Sinv = (3*eye(dsub) - X)/2;
        I1 = (I + I*X)/2;   % I = I * S;
        P1 = (3*P - X*P)/2; % P = Sinv * P;
        % now, P and I are no longer orthogonal, correct
        omega(k) = norm(P1*I1 - eye(d), 'fro');
        if omega(k) > 0.9
            % too much loss of orthogonality, compute explicit inverse
            P1 = (P1*I1)\P1;
        elseif mod(k, 2) == 0
            P1 = 2*P1 - P1*I1*P1;
        else
            I1 = 2*I1 - I1*P1*I1;
        end
        delta(k) = norm(I1 - I, 'fro')/norm(I, 'fro');
        exitFlag = tolerances.test(omega, delta, k);
        k = k + 1;
        I = I1;
        P = P1;
    end
    gen = replab.rep.GenSubRep(rho, gen0.divisionRing, false, I1, P1);
end
