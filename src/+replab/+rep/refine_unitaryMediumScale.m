function Q = refine_unitaryMediumScale(rep, Q, innerIterations, maxIterations)
% Refines an orthogonal basis for subrepresentation of a unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   Q (double(\*,\*)): Orthogonal basis, given as column vectors
%   nInnerIterations (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%
% Returns:
%   double(\*,\*): Refined orthogonal basis
    replab.log(1, 'Unitary refinement: dim(parent) = %d, dim(subrep) = %d', size(Q, 1), size(Q, 2));
    replab.log(1, 'Full commutant projection algorithm');
    replab.log(1, '');
    replab.log(2, ' #iter   dProj    dSpan    ortho');
    replab.log(2, '-----------------------------------');
    iter = 1;
    dsub = size(Q, 2);
    dProj = inf;
    while iter <= maxIterations
        Qprev = Q;
        Ftilde = Q*Q';
        [Foverline err] = rep.commutant.project(Ftilde);
        dProjPrev = dProj;
        dProj = norm(Foverline-Ftilde, 'fro');
        for j = 1:innerIterations
            Q1 = Foverline * Q;
            [Q1, R] = qr(Q1, 0);
            dSpan = norm(Q'*Q1*Q1'*Q - speye(dsub), 'fro')*2;
            replab.log(3, ' (inner)          %6.2E', dSpan);
            Q = Q1;
        end
        ortho = norm(R - speye(dsub), 'fro');
        dSpan = norm(Q'*Qprev*Qprev'*Q - speye(dsub), 'fro')*2;
        replab.log(2, '%6d   %6.2E %6.2E %6.2E', iter, dProj, dSpan, ortho);
        if dProj >= dProjPrev && iter > 3
            break
        end
        iter = iter + 1;
    end
    replab.log(1, 'Stopped after %d iterations with projector error %6.2E', iter, dProj);
end
