function gen1 = unitarize_largeScale(gen, numNonImproving, nSamples, maxIterations)
% Refines a generic non-unitary subrepresentation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to unitarize
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    d = gen.parent.dimension;
    dsub = gen.dimension;
    rep = gen.parent;
    type = [gen.divisionRing '/' rep.field];
    replab.msg(1, 'Unitarization over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_ortho = inf;
    I = gen.injection;
    P = gen.projection;
    ni = 0;
    while iter <= maxIterations
        Iprev = I;
        Pprev = P;
        X = zeros(dsub, dsub);
        for j = 1:nSamples
            g = rep.group.sample;
            switch type
              case {'R/R', 'C/C'}
                rhoI = rep.matrixRowAction(g, Iprev);
              case 'C/R'
                rhoI = rep.matrixRowAction(g, real(Iprev)) + ...
                       rep.matrixRowAction(g, imag(Iprev)) * 1i;
              case 'H/R'
                rhoI = rep.matrixRowAction(g, Iprev.part1) + ...
                       rep.matrixRowAction(g, Iprev.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Iprev.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Iprev.partk) * replab.H.k;
            end
            img = P*rhoI;
            X = X + img*img';
        end
        X = X/trace(X)*dsub;
        % approximation of the matrix square root (see Taylor expansion of sqrt(1+x))
        % S = (eye(dsub) + X)/2;
        % and the square root inverse (see Taylor expansion of 1/sqrt(1+x))
        % Sinv = (3*eye(dsub) - X)/2;
        I = (I + I*X)/2;   % I = I * S;
        P = (3*P - X*P)/2; % P = Sinv * P;
        % now, P and I are no longer orthogonal, correct
        dSpan = norm(P*I - eye(dsub), 'fro');
        if dSpan > 0.9
            % too much loss of orthogonality, compute explicit inverse
            P = (P*I)\P;
        elseif mod(iter, 2) == 0
            P = 2*P - P*I*P;
        else
            I = 2*I - I*P*I;
        end
        ortho = norm(X - eye(dsub), 'fro');
        if ortho >= min_ortho && ortho < 1
            ni = ni + 1;
            replab.msg(2, '%6d   %6.2E %6.2E (#%d non improving)', iter, dSpan, ortho, ni);
            if ni > numNonImproving
                break
            end
        else
            min_ortho = ortho;
            replab.msg(2, '%6d   %6.2E %6.2E', iter, dSpan, ortho);
        end
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    gen1 = replab.rep.GenSubRep(rep, gen.divisionRing, false, I, P);
end
