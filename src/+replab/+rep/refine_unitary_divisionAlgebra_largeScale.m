function basis = refine_unitary_divisionAlgebra_largeScale(rep, basis0, divisionAlgebraName, numNonImproving, nSamples, maxIterations, Qo)
% Refines an orthogonal basis for subrepresentation of a unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation, must be unitary
%   basis0 (double(\*,\*)): Orthogonal basis, given as column vectors
%   divisionAlgebraName ('', 'complex', 'quaternion.rep'): Division algebra to preserve
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Qo (double(D,do)): Basis matrix prescribing orthogonality, used in isotypic components
%
% Returns:
%   double(\*,\*): Refined orthogonal basis
    d = size(basis0, 1);
    dsub = size(basis0, 2);
    replab.msg(1, 'Unitary refinement: dim(parent) = %d, dim(subrep) = %d', d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
    switch divisionAlgebraName
      case ''
        Q0 = basis0;
        dda = 1;
      case 'complex'
        A0 = basis0(:,1:2:end)/sqrt(2);
        B0 = basis0(:,2:2:end)/sqrt(2);
        Q0 = A0 + 1i*B0;
        dda = 2;
      case 'quaternion.rep'
        i = replab.H.i;
        j = replab.H.j;
        k = replab.H.k;
        U = [1+i+j+k; 1-i+j-k; 1-i-j+k; 1+i-j-k]/4;
        Q0 = basis0 * kron(eye(dsub/4), U);
        dda = 4;
    end
    Q = Q0;
    ni = 0;
    while iter <= maxIterations
        Q1 = zeros(size(Q));
        for j = 1:nSamples
            g = rep.group.sample;
            switch divisionAlgebraName
              case ''
                Qrho = rep.matrixColAction(g, Q');
                rhoQ = rep.matrixRowAction(g, Q);
              case 'complex'
                Qrho = rep.matrixColAction(g, real(Q)') + ...
                       rep.matrixColAction(g, -imag(Q)') * 1i;
                rhoQ = rep.matrixRowAction(g, real(Q)) + ...
                       rep.matrixRowAction(g, imag(Q)) * 1i;
              case 'quaternion.rep'
                Qrho = rep.matrixColAction(g, Q.part1') + ...
                       rep.matrixColAction(g, -Q.parti') * replab.H.i + ...
                       rep.matrixColAction(g, -Q.partj') * replab.H.j + ...
                       rep.matrixColAction(g, -Q.partk') * replab.H.k;
                rhoQ = rep.matrixRowAction(g, Q.part1) + ...
                       rep.matrixRowAction(g, Q.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Q.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Q.partk) * replab.H.k;
            end
            Q1 = Q1 + rhoQ * (Qrho * Q);
        end
        Q1 = Q1/nSamples;
        if ~isempty(Qo)
            Q1 = Q1 - Qo * (Qo' * Q1);
        end
        [Q1, R] = qr(Q1, 0);
        ortho = norm(R - speye(dsub/dda), 'fro');
        dSpan = norm(Q1'*Q*Q'*Q1 - speye(dsub/dda), 'fro');
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
        Q = Q1;
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    switch divisionAlgebraName
      case ''
        basis = Q;
      case 'complex'
        A = real(Q);
        B = imag(Q);
        basis = zeros(d, dsub);
        basis(:,1:2:end) = A*sqrt(2);
        basis(:,2:2:end) = B*sqrt(2);
      case 'quaternion.rep'
        A = Q.part1;
        B = Q.parti;
        C = Q.partj;
        D = Q.partk;
        basis = zeros(d, dsub);
        basis(:,1:4:end) = (A+B+C+D);
        basis(:,2:4:end) = (A-B-C+D);
        basis(:,3:4:end) = (A+B-C-D);
        basis(:,4:4:end) = (A-B+C-D);
    end
end
