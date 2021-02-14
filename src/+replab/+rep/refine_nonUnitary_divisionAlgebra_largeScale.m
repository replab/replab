function [inj, prj, I, P] = refine_nonUnitary_divisionAlgebra_largeScale(rep, inj0, prj0, divisionAlgebraName, numNonImproving, nSamples, maxIterations, Ibo, Pbo)
% Refines an injection/projection pair for subrepresentation of a possibly non-unitary representation
%
% Args:
%   rep (`+replab.Rep`): Parent representation
%   inj0 (double(\*,\*)): Injection map matrix
%   prj0 (double(\*,\*)): Projection map matrix
%   divisionAlgebraName ('', 'complex', 'quaternion.rep'): Division algebra to preserve
%   numNonImproving (integer): See `+replab.SubRep.refine`
%   nSamples (integer): See `+replab.SubRep.refine`
%   maxIterations (integer): See `+replab.SubRep.refine`
%   Ibo (double(D,do)): Injection map matrix prescribing biorthogonality
%   Pbo (double(do,D)): Projection map matrix prescribing biorthogonality
%
% Returns
% -------
%   inj: double(\*,\*)
%     Refined injection map
%   prj: double(\*,\*)
%     Refined projection map
    d = rep.dimension;
    dsub = size(inj0, 2);
    replab.msg(1, 'Nonunitary refinement: dim(parent) = %d, dim(subrep) = %d', d, dsub);
    replab.msg(1, 'Large-scale algorithm with %d samples/iteration', nSamples);
    replab.msg(1, '');
    replab.msg(2, ' #iter   dSpan    ortho');
    replab.msg(2, '--------------------------');
    iter = 1;
    min_dSpan = inf;
    switch divisionAlgebraName
      case ''
        dda = 1;
        I0 = inj0;
        P0 = prj0;
      case 'complex'
        % the complex subrepresentation would be defined by
        % I = A+1i*B
        % P = E-1i*F
        % where A,B is interleaved in I, and E,F is interleaved in P
        dda = 2;
        A0 = inj0(:,1:2:end)/sqrt(2);
        B0 = inj0(:,2:2:end)/sqrt(2);
        E0 = prj0(1:2:end,:)/sqrt(2);
        F0 = prj0(2:2:end,:)/sqrt(2);
        I0 = (A0 + 1i*B0);
        P0 = (E0 - 1i*F0);
      case 'quaternion.rep'
        dda = 4;
        i = replab.H.i;
        j = replab.H.j;
        k = replab.H.k;
        U = [1+i+j+k; 1-i+j-k; 1-i-j+k; 1+i-j-k]/4;
        I0 = inj0 * kron(eye(dsub/4), U);
        P0 = kron(eye(dsub/4), U)' * prj0;
    end
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
            switch divisionAlgebraName
              case ''
                Prho = rep.matrixColAction(g, Pprev);
                rhoI = rep.matrixRowAction(g, Iprev);
              case 'complex'
                Prho = rep.matrixColAction(g, real(Pprev)) + ...
                       rep.matrixColAction(g, imag(Pprev)) * 1i;
                rhoI = rep.matrixRowAction(g, real(Iprev)) + ...
                       rep.matrixRowAction(g, imag(Iprev)) * 1i;
              case 'quaternion.rep'
                Prho = rep.matrixColAction(g, Pprev.part1) + ...
                       rep.matrixColAction(g, Pprev.parti) * replab.H.i + ...
                       rep.matrixColAction(g, Pprev.partj) * replab.H.j + ...
                       rep.matrixColAction(g, Pprev.partk) * replab.H.k;
                rhoI = rep.matrixRowAction(g, Iprev.part1) + ...
                       rep.matrixRowAction(g, Iprev.parti) * replab.H.i + ...
                       rep.matrixRowAction(g, Iprev.partj) * replab.H.j + ...
                       rep.matrixRowAction(g, Iprev.partk) * replab.H.k;
            end
            I = I + rhoI * (Prho * Iprev);
            P = P + (Pprev * rhoI) * Prho;
        end
        if ~isempty(Ibo) && ~isempty(Pbo)
            I = I - Ibo * (Pbo * I);
            P = P - (P * Ibo) * Pbo;
        end
        I = I / (P0 * I);
        P = (P * I) \ P;
        [U, ~] = qr(I, 0);
        dSpan = norm(U'*Uprev*Uprev'*U - speye(dsub/dda), 'fro');
        ortho = norm(P * Iprev - speye(dsub/dda), 'fro');
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
        iter = iter + 1;
    end
    replab.msg(1, 'Stopped after %d iterations with span delta %6.2E', iter, dSpan);
    switch divisionAlgebraName
      case ''
        inj = I;
        prj = P;
      case 'complex'
        A = real(I);
        B = imag(I);
        E = real(P);
        F = -imag(P);
        inj = zeros(d, dsub);
        prj = zeros(dsub, d);
        inj(:,1:2:end) = A*sqrt(2);
        inj(:,2:2:end) = B*sqrt(2);
        prj(1:2:end,:) = E*sqrt(2);
        prj(2:2:end,:) = F*sqrt(2);
        % final correction: why?
        prj = (prj * inj) \ prj;
      case 'quaternion.rep'
        A = I.part1;
        B = I.parti;
        C = I.partj;
        D = I.partk;
        E = P.part1;
        F = P.parti;
        G = P.partj;
        H = P.partk;
        inj = zeros(d, dsub);
        prj = zeros(dsub, d);
        inj(:,1:4:end) = (A+B+C+D);
        inj(:,2:4:end) = (A-B-C+D);
        inj(:,3:4:end) = (A+B-C-D);
        inj(:,4:4:end) = (A-B+C-D);
        prj(1:4:end,:) = (E-F-G-H);
        prj(2:4:end,:) = (E+F+G-H);
        prj(3:4:end,:) = (E-F+G+H);
        prj(4:4:end,:) = (E+F-G+H);
        % final correction: why?
        prj = (prj * inj) \ prj;
    end
end
