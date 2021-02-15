function gen1 = harmonize_unitary_largeScale(gen, gen0, numNonImproving, nSamples, maxIterations, Qo)
% Refines a generic non-unitary subrepresentation
%
% Args:
%   gen (`+replab.GenSubRep`): Generic subrepresentation to harmonize
%   gen0 (`+replab.GenSubRep`): Reference generic subrepresentation
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Qo (double(D,do)): Injection map matrix prescribing orthogonality
%
% Returns:
%   `+replab.GenSubRep`: Refined generic subrepresentation
    d = gen.parent.dimension;
    dsub = gen.dimension;
    rep = gen.parent;
    type = [gen.divisionRing '/' rep.field];
    replab.msg(1, 'Unitary harmonization over %s: dim(parent) = %d, dim(subrep) = %d', gen.divisionRing, d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_ortho = inf;
    Q0 = gen0.injection;
    Q = gen.injection;
    assert(rep.isUnitary);
    if ~gen.mapsAreAdjoint
        [Q, ~] = replab.numerical.qr(Q);
    end
    ni = 0;
    while iter <= maxIterations
        Q1 = zeros(size(Q));
        for j = 1:nSamples
            g = rep.group.sample;
            switch type
              case {'R/R', 'C/C'}
                rhoQ0 = rep.matrixRowAction(g, Q0);
                rhoQ = rep.matrixRowAction(g, Q);
              case 'C/R'
                rhoQ0 = rep.matrixRowAction(g, real(Q0)) + ...
                        rep.matrixRowAction(g, imag(Q0)) * 1i;
                rhoQ = rep.matrixRowAction(g, real(Q)) + ...
                       rep.matrixRowAction(g, imag(Q)) * 1i;
              case 'H/R'
                rhoQ0 = rep.matrixRowAction(g, Q0.part1) + ...
                        rep.matrixRowAction(g, Q0.parti) * replab.H.i + ...
                        rep.matrixRowAction(g, Q0.partj) * replab.H.j + ...
                        rep.matrixRowAction(g, Q0.partk) * replab.H.k;
                rhoQ = rep.matrixRowAction(g, Q.part1) + ...
                       rep.matrixRowAction(g, Q.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Q.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Q.partk) * replab.H.k;
            end
            Q1 = Q1 + rhoQ * (rhoQ0' * Q0);
        end
        if ~isempty(Qo)
            Q1 = Q1 - Qo * (Qo' * Q1);
        end
        f = trace(Q1*Q1')/dsub;
        Q1 = Q1/sqrt(f);
        dSpan = norm(Q1'*Q - eye(dsub), 'fro');
        ortho = norm(Q1'*Q1 - eye(dsub), 'fro');
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
        Q = Q1;
        iter = iter + 1;
    end
    % One step of Newton iteration to force unitary (helps!)
    N = Q'*Q;
    P = Q*N/2;
    Q = 2*Q + P*N - 3*P;
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    gen1 = replab.rep.GenSubRep(rep, gen.divisionRing, true, Q, Q');
end
