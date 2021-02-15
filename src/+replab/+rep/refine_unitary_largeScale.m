function gen1 = refine_unitary_largeScale(gen, numNonImproving, nSamples, maxIterations, Qo)
% Refines a generic unitary subrepresentation
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
    replab.msg(1, 'Unitary refinement over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
    rep = gen.parent;
    Q0 = gen.injection;
    if ~gen.mapsAreAdjoint
        [Q0, ~] = replab.numerical.qr(Q0);
    end
    Q = Q0;
    ni = 0;
    while iter <= maxIterations
        Q1 = zeros(size(Q));
        for j = 1:nSamples
            g = rep.group.sample;
            switch gen.divisionRing
              case 'R'
                Qrho = rep.matrixColAction(g, Q');
                rhoQ = rep.matrixRowAction(g, Q);
              case 'C'
                Qrho = rep.matrixColAction(g, real(Q)') + ...
                       rep.matrixColAction(g, -imag(Q)') * 1i;
                rhoQ = rep.matrixRowAction(g, real(Q)) + ...
                       rep.matrixRowAction(g, imag(Q)) * 1i;
              case 'H'
                Qrho = replab.H(rep.matrixColAction(g, Q.part1'), ...
                                rep.matrixColAction(g, -Q.parti'), ...
                                rep.matrixColAction(g, -Q.partj'), ...
                                rep.matrixColAction(g, -Q.partk'));
                rhoQ = replab.H(rep.matrixRowAction(g, Q.part1), ...
                                rep.matrixRowAction(g, Q.parti), ...
                                rep.matrixRowAction(g, Q.partj), ...
                                rep.matrixRowAction(g, Q.partk));
            end
            Q1 = Q1 + rhoQ * (Qrho * Q);
        end
        Q1 = Q1/nSamples;
        if ~isempty(Qo)
            Q1 = Q1 - Qo * (Qo' * Q1);
        end
        Qbefore = Q1;
        [Q1, ~] = replab.numerical.qr(Q1);
        %nIters = 1;
        %for j = 1:nIters
        %    N = Q1'*Q1;
        %    P = Q1*N/2;
        %    Q1 = 2*Q1 + P*N - 3*P;
        %end
        %Q2'*Q1
        ortho = norm(Qbefore'*Q1 - speye(dsub), 'fro');
        dSpan = norm(Q1'*Q*Q'*Q1 - speye(dsub), 'fro');
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
        Q = Q1;
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    gen1 = replab.rep.GenSubRep(rep, gen.divisionRing, true, Q, Q');
end
