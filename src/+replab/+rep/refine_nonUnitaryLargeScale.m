function [I, P] = refine_nonUnitaryLargeScale(rep, I0, P0, numNonImproving, nSamples, maxIterations, Ibo, Pbo)
% Refines an injection/projection pair for subrepresentation of a possibly non-unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   I0 (double(\*,\*)): Injection map matrix
%   P0 (double(\*,\*)): Projection map matrix
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Ibo (double(D,do)): Injection map matrix prescribing biorthogonality
%   Pbo (double(do,D)): Projection map matrix prescribing biorthogonality
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
    replab.log(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.log(1, '');
    replab.log(2, ' #iter   dSpan    ortho');
    replab.log(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
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
            I = I + rep.matrixRowAction(g, Iprev) * (rep.matrixColAction(g, Pprev) * Iprev);
            P = P + (Pprev * rep.matrixRowAction(g, Iprev)) * rep.matrixColAction(g, Pprev);
        end
        if ~isempty(Ibo) && ~isempty(Pbo)
            I = I - Ibo * (Pbo * I);
            P = P - (P * Ibo) * Pbo;
        end
        I = I / (P0 * I);
        P = (P * I) \ P;
        [U, ~] = qr(I, 0);
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub), 'fro')*2;
        ortho = norm(P * Iprev - speye(dsub), 'fro');
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
        iter = iter + 1;
    end
    replab.log(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
end
