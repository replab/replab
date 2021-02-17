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
    replab.msg(2, ' #iter   ortho    delta');
    replab.msg(2, '-----------------------');
    iter = 1;
    P = I'; % as good as a guess as anything else
    [P, min_ortho] = replab.rep.biorthoStepP(I, P);
    ni = 0;
    while iter <= maxIterations
        Pprev = P;
        P = zeros(size(Pprev));
        for j = 1:nSamples
            g = rep.group.sample;
            P = P + (Pprev * rep.matrixRowAction(g, I)) * rep.matrixColAction(g, Pprev);
        end
        f = trace(P*I)/dsub;
        P = P/f;
        delta = norm(P - Pprev, 'fro');
        [P, ortho] = replab.rep.biorthoStepP(I, P);
        if ortho > min_ortho
            ni = ni + 1;
            replab.msg(2, '%6d   %6.2E %6.2E (#%d non improving)', iter, ortho, delta, ni);
            if ni > numNonImproving
                break
            end
        else
            min_ortho = ortho;
            replab.msg(2, '%6d   %6.2E %6.2E', iter, ortho, delta);
        end
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with ortho %6.2E', iter, ortho);
end
