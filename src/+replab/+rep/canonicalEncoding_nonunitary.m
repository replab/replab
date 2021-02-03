function subs = canonicalEncoding_nonunitary(sub, iterator)
% Identifies the real representation(s) present in a real subrepresentation encoding a pair of conjugate irreps
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation with `+replab.SubRep.encodesIrreducibleComplexPair` set to true
%   iterator (`+replab.+domain.SamplesIterator`): Iterator in the sequence of parent commutant samples
%
% Returns:
%   cell(1,\*) of `+replab.SubRep`: Real irreps with `+replab.Rep.frobeniusSchurIndicator` computed and `+replab.Rep.isDivisionAlgebraCanonical` set to true
    assert(sub.overR);
    assert(sub.encodesIrreducibleComplexPair);
    d = sub.dimension;
    I = sub.injection('double/sparse');
    P = sub.projection('double/sparse');
    d = sub.dimension;
    A = I(:, 1:2:d);
    B = I(:, 2:2:d);
    C = P(1:2:d, :);
    D = P(2:2:d, :);
    U = [A+1i*B A-1i*B];
    V = [C-1i*D;C+1i*D];
    S = iterator.next;
    X = V*S*U;
    J = X(1:d/2, d/2+1:d);
    cJ = X(d/2+1:d, 1:d/2);
    cJJ = conj(J)*J;
    lambda = trace(cJJ)/(d/2);
    err = cJJ - lambda*eye(d/2);
    tol = replab.globals.doubleEigTol;
    if abs(lambda) < tol
        % complex-type representation
        fs = 0;
        sub.cache('frobeniusSchurIndicator', 0, '==');
        sub.cache('isDivisionAlgebraCanonical', true, '==');
        subs = {sub};
    elseif lambda > 0
        % pair of real-type representations
    else % lambda < 0
         % quaternion-type representation
        assert(mod(d, 4) == 0);
        [U1, D1] = eig(J, 'vector');
        [~, I1] = sort(real(D1));
        I1 = I1(:)';
        assert(all(D1(I1) + D1(fliplr(I1)) < tol));
        I1 = [I1(1:d/4) fliplr(I1(d/4+1:end))];
        D1 = D1(I1);
        U1 = U1(:,I1);
        V1 = inv(U1);
        U2 = U1(:,1:d/4);
        V2 = V1(1:d/4,:);
        UU = blkdiag([U2 conj(U2)], [U2 conj(U2)]);
        VV = blkdiag([V2; conj(V2)], [V2; conj(V2)]);
        XX=VV*X*UU;
        v = diag(XX(1:d/2,d/2+1:end));
        v = v(1:d/4);
        v = v/abs(v(1));
        w = sqrt(v);
        w = w(:).';
        %w = [w 1i*w];
        W = diag([w 1i*w conj(w) -1i*conj(w)]);
        s1 = [1 0; 0 1];
        sy = [0 -1i; 1i 0];
        P = [1:2:d/2 2:2:d/2 d/2+(1:2:d/2) d/2+(2:2:d/2)];
        P = sparse(P(:)', 1:d, ones(1, d), d, d);
        W1 = kron([s1-sy s1+sy; -s1-sy s1-sy]/2, eye(d/4));
        XX = W1*P*conj(W) * VV*V*S *U*UU * W*P'*W1';
        XX1 = W1*P*conj(W) * VV*V*iterator.next *U*UU * W*P'*W1';
        %XX1 = diag([conj(w) w]) * VV*V* iterator.next *U*UU * diag([w conj(w)]);
        % has the form [0 L; conj(L) 0] with L = blkdiag(L1, -L1)
        error('asd');
% $$$         Y = V2*V*S*U*U2;
% $$$         D2 = ones(1, d);
% $$$         for i = 2:d/2
% $$$             f = Y(1,d/2+1)/Y(i,d/2+i);
% $$$             D2(d/2+i) = f;
% $$$         end
% $$$         Y = diag(1./D2)*V2*V*S*U*U2*diag(D2);
% $$$         f = sqrt(Y(d/2+1,1)/Y(1,d/2+1));
% $$$         D2(d/2+1:end) = f * D2(d/2+1:end);
% $$$         Y = diag(1./D2)*V2*V*S*U*U2*diag(D2);
        sz = [1 0; 0 -1];
        sx = [0 1; 1 0];
        W11 = sz+1i*sx;
        W12 = sz-1i*sx;
        W21 = -1i*sz-sx;
        W22 = -1i*sz+sx;
        I = speye(d/4);
        W = [kron(I, W11) kron(I, W12)
             kron(I, W21) kron(I, W22)];
% $$$         Y = W*diag(1./D2)*V2*V*S*U*U2*diag(D2)*W'
    end
end