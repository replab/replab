function gen1 = refine_unitary_mediumScale(gen, innerIterations, maxIterations, Qo)
% Refines a generic unitary representation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to refine
%   nInnerIterations (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Qo (double(D,do)): Basis matrix prescribing orthogonality, used in isotypic components
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    d = gen.parent.dimension;
    dsub = gen.dimension;
    rep = gen.parent;
    replab.msg(1, 'Unitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Full commutant projection algorithm');
    replab.msg(1, '');
    replab.msg(2, ' #iter   dProj    dSpan    ortho');
    replab.msg(2, '-----------------------------------');
    iter = 1;
    Q0 = gen.injection;
    if ~gen.mapsAreAdjoint
        [Q0, ~] = replab.numerical.qr(Q0);
    end
    Q = Q0;
    dProj = inf;
    while iter <= maxIterations
        Qprev = Q;
        Ftilde = Q*Q';
        switch gen.divisionRing
          case 'R'
            [Foverline, err] = rep.commutant.project(Ftilde);
          case 'C'
            [F1, e1] = rep.commutant.project(real(Ftilde));
            [Fi, ei] = rep.commutant.project(imag(Ftilde));
            Foverline = F1 + Fi*1i;
            err = sqrt(e1^2 + ei^2);
          case 'H'
            [F1, e1] = rep.commutant.project(Ftilde.part1);
            [Fi, ei] = rep.commutant.project(Ftilde.parti);
            [Fj, ej] = rep.commutant.project(Ftilde.partj);
            [Fk, ek] = rep.commutant.project(Ftilde.partk);
            Foverline = replab.H(F1, Fi, Fj, Fk);
            err = sqrt(e1^2 + ei^2 + ej^2 + ek^2);
        end
        dProjPrev = dProj;
        dProj = norm(Foverline-Ftilde, 'fro');
        for j = 1:innerIterations
            Q1 = Foverline * Q;
            if ~isempty(Qo)
                Q1 = Q1 - Qo * (Qo' * Q1);
            end
            [Q1, R] = replab.numerical.qr(Q1);
            dSpan = norm(Q'*Q1*Q1'*Q - speye(dsub), 'fro')*2;
            replab.msg(3, ' (inner)          %6.2E', dSpan);
            Q = Q1;
        end
        ortho = norm(R - speye(dsub), 'fro');
        dSpan = norm(Q'*Qprev*Qprev'*Q - speye(dsub), 'fro')*2;
        replab.msg(2, '%6d   %6.2E %6.2E %6.2E', iter, dProj, dSpan, ortho);
        if dProj >= dProjPrev && iter > 3
            break
        end
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with projector error %6.2E', iter, dProj);
    gen1 = replab.rep.GenSubRep(rep, gen.divisionRing, true, Q, Q');
end
