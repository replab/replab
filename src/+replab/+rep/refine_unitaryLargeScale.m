function Q = refine_unitaryLargeScale(rep, Q, numNonImproving, nSamples, maxIterations)
% Refines an orthogonal basis for subrepresentation of a unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   Q (double(\*,\*)): Orthogonal basis, given as column vectors
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%
% Returns:
%   double(\*,\*): Refined orthogonal basis
    replab.log(1, 'Unitary refinement: dim(parent) = %d, dim(subrep) = %d', size(Q, 1), size(Q, 2));
    replab.log(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.log(1, '');
    replab.log(2, ' #iter   dSpan    ortho');
    replab.log(2, '--------------------------');
    iter = 1;
    dsub = size(Q, 2);
    min_dSpan = inf;
    ni = 0;
    while iter <= maxIterations
        Q1 = zeros(size(Q));
        for j = 1:nSamples
            g = rep.group.sample;
            Q1 = Q1 + rep.matrixRowAction(g, Q) * (rep.matrixColAction(g, Q') * Q);
        end
        Q1 = Q1/nSamples;
        [Q1, R] = qr(Q1, 0);
        ortho = norm(R - speye(dsub), 'fro');
        dSpan = norm(Q1'*Q*Q'*Q1 - speye(dsub), 'fro')*2;
        if dSpan > min_dSpan
            ni = ni + 1;
            replab.log(2, '%6d   %6.2E %6.2E (#%d non improving)', iter, dSpan, ortho, ni);
            if ni > numNonImproving
                break
            end
        else
            min_dSpan = dSpan;
            replab.log(2, '%6d   %6.2E %6.2E', iter, dSpan, ortho);
        end
        Q = Q1;
    end
    replab.log(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
end
