function P = findProjection_largeScale(rep, I, numNonImproving, nSamples, maxIterations)
% Finds a projection map for a subrepresentation defined by an injection map
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   I (double(\*,\*)): Injection map matrix
%   numNonImproving (integer): Number of non-improving steps before stopping the large-scale algorithm
%   nSamples (integer): Number of samples to use in the large-scale version of the algorithm
%   maxIterations (integer): Maximum number of (outer) iterations, default ``1000``
%
% Returns:
%   double(\*,\*): Projection map
    d = rep.dimension;
    dsub = size(I, 2);
    replab.msg(1, 'Projection map search: dim(parent) = %d, dim(subrep) = %d', d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
    P = I';
    ni = 0;
    [U, ~] = replab.numerical.qr(P');
    while iter <= maxIterations
        Pprev = P;
        Uprev = U;
        P = zeros(size(Pprev));
        for j = 1:nSamples
            g = rep.group.sample;
            P = P + (Pprev * rep.matrixRowAction(g, I)) * rep.matrixColAction(g, Pprev);
        end
        P = (P * I) \ P;
        [U, ~] = replab.numerical.qr(P');
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub), 'fro')*2;
        ortho = norm(P * I - speye(dsub), 'fro');
        if dSpan > min_dSpan
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
end
