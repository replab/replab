function gen1 = refine_nonUnitary_genSubRep_largeScale(gen, numNonImproving, nSamples, maxIterations, Ibo, Pbo)
% Refines a generic non-unitary subrepresentation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to refine
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
    replab.msg(1, 'Nonunitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
    rep = gen.parent;
    I0 = gen.injection;
    P0 = gen.projection;
    I = I0;
    P = P0;
    ni = 0;
    [U, ~] = qr(I, 0);
    while iter <= maxIterations
        Iprev = I;
        Pprev = P;
        Uprev = U;
        I = zeros(size(Iprev));
        P = zeros(size(Pprev));
        for j = 1:nSamples
            g = rep.group.sample;
            switch gen.divisionRing
              case 'R'
                Prho = rep.matrixColAction(g, Pprev);
                rhoI = rep.matrixRowAction(g, Iprev);
              case 'C'
                Prho = rep.matrixColAction(g, real(Pprev)) + ...
                       rep.matrixColAction(g, imag(Pprev)) * 1i;
                rhoI = rep.matrixRowAction(g, real(Iprev)) + ...
                       rep.matrixRowAction(g, imag(Iprev)) * 1i;
              case 'H'
                Prho = rep.matrixColAction(g, Pprev.part1) + ...
                       rep.matrixColAction(g, Pprev.parti) * replab.H.i + ...
                       rep.matrixColAction(g, Pprev.partj) * replab.H.j + ...
                       rep.matrixColAction(g, Pprev.partk) * replab.H.k;
                rhoI = rep.matrixRowAction(g, Iprev.part1) + ...
                       rep.matrixRowAction(g, Iprev.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Iprev.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Iprev.partk) * replab.H.k;
            end
            I = I + rhoI * (Prho * Iprev);
            P = P + (Pprev * rhoI) * Prho;
        end
        if ~isempty(Ibo) && ~isempty(Pbo)
            I = I - Ibo * (Pbo * I);
            P = P - (P * Ibo) * Pbo;
        end
        I = I / (P0 * I);
        P = (P * I) \ P;
        [U, ~] = qr(I, 0);
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub), 'fro');
        ortho = norm(P * Iprev - speye(dsub), 'fro');
        if dSpan >= min_dSpan
            ni = ni + 1;
            replab.msg(2, '%6d   %6.2E %6.2E (#%d non improving)', iter, dSpan, ortho, ni);
            if ni > numNonImproving
                break
            end
        else
            min_dSpan = dSpan;
            replab.msg(2, '%6d   %6.2E %6.2E', iter, dSpan, ortho);
        end
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    gen1 = replab.GenSubRep(rep, gen.divisionRing, false, I, P);
end
