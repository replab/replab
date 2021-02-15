function gen1 = refine_nonUnitary_mediumScale(gen, innerIterations, maxIterations, Ibo, Pbo)
% Refines a generic non-unitary representation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to refine
%   nInnerIterations (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Ibo (double(D,do)): Injection map matrix prescribing biorthogonality (in the parent representation field)
%   Pbo (double(do,D)): Projection map matrix prescribing biorthogonality (in the parent representation field)
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    d = gen.parent.dimension;
    dsub = gen.dimension;
    rep = gen.parent;
    replab.msg(1, 'Nonunitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Using the full commutant projection algorithm');
    replab.msg(1, '');
    replab.msg(2, ' #iter   dProj    dSpan    ortho');
    replab.msg(2, '-----------------------------------');
    iter = 1;
    dProj = inf;
    I0 = gen.injection;
    P0 = gen.projection;
    I = I0;
    P = P0;
    [U, ~] = replab.numerical.qr(I);
    while iter <= maxIterations
        Iprev = I;
        Pprev = P;
        Uprev = U;
        Ftilde = I*P;
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
            I1 = Foverline * I;
            P1 = P * Foverline;
            if ~isempty(Ibo) && ~isempty(Pbo)
                I1 = I1 - Ibo * (Pbo * I1);
                P1 = P1 - (P1 * Ibo) * Pbo;
            end
            I1 = I1 / (P0 * I1);
            P1 = (P1 * I1) \ P1;
            [U1, ~] = replab.numerical.qr(I1);
            ortho = norm(P1 * I0 - speye(dsub), 'fro');
            dSpan = norm(U'*U1*U1'*U - speye(dsub), 'fro')*2;
            replab.msg(3, ' (inner)          %6.2E %6.2E', dSpan, ortho);
            U = U1;
            I = I1;
            P = P1;
        end
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub), 'fro')*2;
        replab.msg(2, '%6d   %6.2E %6.2E %6.2E', iter, dProj, dSpan, ortho);
        if dProj >= dProjPrev && iter > 3
            break
        end
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with projector error %6.2E', iter, dProj);
    gen1 = replab.rep.GenSubRep(rep, gen.divisionRing, false, I, P);
end
