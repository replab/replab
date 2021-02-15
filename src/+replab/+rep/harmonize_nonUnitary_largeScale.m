function gen1 = harmonize_nonUnitary_largeScale(gen, gen0, numNonImproving, nSamples, maxIterations, Ibo, Pbo)
% Refines a generic non-unitary subrepresentation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to harmonize
%   gen0 (`+replab.GenSubRep`): Reference generic subrepresentation
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Ibo (double(D,do)): Injection map matrix prescribing biorthogonality
%   Pbo (double(do,D)): Projection map matrix prescribing biorthogonality
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    d = gen.parent.dimension;
    dsub = gen.dimension;
    rep = gen.parent;
    type = [gen.divisionRing '/' rep.field];
    replab.msg(1, 'Nonunitary harmonization over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_ortho = inf;
    I0 = gen0.injection;
    P0 = gen0.projection;
    I = gen.injection;
    P = gen.projection;
    ni = 0;
    while iter <= maxIterations
        Iprev = I;
        Pprev = P;
        I = zeros(size(Iprev));
        for j = 1:nSamples
            g = rep.group.sample;
            switch type
              case {'R/R', 'C/C'}
                P0rho = rep.matrixColAction(g, P0);
                rhoI = rep.matrixRowAction(g, Iprev);
              case 'C/R'
                P0rho = rep.matrixColAction(g, real(P0)) + ...
                        rep.matrixColAction(g, imag(P0)) * 1i;
                rhoI = rep.matrixRowAction(g, real(Iprev)) + ...
                       rep.matrixRowAction(g, imag(Iprev)) * 1i;
              case 'H/R'
                P0rho = rep.matrixColAction(g, P0.part1) + ...
                        rep.matrixColAction(g, P0.parti) * replab.H.i + ...
                        rep.matrixColAction(g, P0.partj) * replab.H.j + ...
                        rep.matrixColAction(g, P0.partk) * replab.H.k;
                rhoI = rep.matrixRowAction(g, Iprev.part1) + ...
                       rep.matrixRowAction(g, Iprev.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Iprev.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Iprev.partk) * replab.H.k;
            end
            I = I + rhoI * (P0rho * I0);
        end
        P = zeros(size(Pprev));
        for j = 1:nSamples
            g = rep.group.sample;
            switch type
              case {'R/R', 'C/C'}
                rhoI0 = rep.matrixRowAction(g, I0);
                Prho = rep.matrixColAction(g, Pprev);
              case 'C/R'
                rhoI0 = rep.matrixRowAction(g, real(I0)) + ...
                        rep.matrixRowAction(g, imag(I0)) * 1i;
                Prho = rep.matrixColAction(g, real(Pprev)) + ...
                       rep.matrixColAction(g, imag(Pprev)) * 1i;
              case 'H/R'
                rhoI0 = rep.matrixRowAction(g, I0.part1) + ...
                        rep.matrixRowAction(g, I0.parti) * replab.H.i + ...
                        rep.matrixRowAction(g, I0.partj) * replab.H.j + ...
                        rep.matrixRowAction(g, I0.partk) * replab.H.k;
                Prho = rep.matrixColAction(g, Pprev.part1) + ...
                       rep.matrixColAction(g, Pprev.parti) * replab.H.i + ...
                       rep.matrixColAction(g, Pprev.partj) * replab.H.j + ...
                       rep.matrixColAction(g, Pprev.partk) * replab.H.k;
            end
            P = P + (P0 * rhoI0) * Prho;
        end
        if ~isempty(Ibo) && ~isempty(Pbo)
            I = I - Ibo * (Pbo * I);
            P = P - (P * Ibo) * Pbo;
        end
        f = trace(P*I)/dsub;
        sf = sqrt(f);
        I = I/sf;
        P = P/sf;
        dSpan = norm(P * Iprev - speye(dsub), 'fro');
        ortho = norm(P*I - eye(dsub), 'fro');
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
