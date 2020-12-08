function [I, P] = refine_nonUnitaryMediumScale(rep, I0, P0, innerIterations, maxIterations)
% Refines an injection/projection pair for subrepresentation of a possibly non-unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   I0 (double(\*,\*)): Injection map matrix
%   P0 (double(\*,\*)): Projection map matrix
%   nInnerIterations (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%
% Returns
% -------
%   I: double(\*,\*)
%     Refined injection map
%   P: double(\*,\*)
%     Refined projection map
    d = rep.dimension;
    dsub = size(I0, 2);
    replab.log(1, 'Nonunitary refinement: dim(parent) = %d, dim(subrep) = %d', d, dsub);
    replab.log(1, 'Using the full commutant projection algorithm');
    replab.log(1, '');
    replab.log(2, ' #iter   dProj    dSpan    ortho');
    replab.log(2, '-----------------------------------');
    iter = 1;
    dProj = inf;
    I = I0;
    P = P0;
    [U, ~] = qr(I, 0);
    while iter <= maxIterations
        Iprev = I;
        Pprev = P;
        Uprev = U;
        Ftilde = I*P;
        [Foverline err] = rep.commutant.project(Ftilde);
        dProjPrev = dProj;
        dProj = norm(Foverline-Ftilde, 'fro');
        for j = 1:innerIterations
            I1 = Foverline * I;
            P1 = P * Foverline;
            I1 = I1 / (P0 * I1);
            P1 = (P1 * I1) \ P1;
            [U1, ~] = qr(I1, 0);
            ortho = norm(P1 * I0 - speye(dsub), 'fro');
            dSpan = norm(U'*U1*U1'*U - speye(dsub), 'fro')*2;
            replab.log(3, ' (inner)          %6.2E %6.2E', dSpan, ortho);
            U = U1;
            I = I1;
            P = P1;
        end
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub), 'fro')*2;
        replab.log(2, '%6d   %6.2E %6.2E %6.2E', iter, dProj, dSpan, ortho);
        if dProj >= dProjPrev && iter > 3
            break
        end
        iter = iter + 1;
    end
    replab.log(1, 'Stopped after %d iterations with projector error %6.2E', iter, dProj);
end
