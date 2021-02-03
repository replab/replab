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
        f = trace(conj(J)*J)/(d/2);
        J = J / sqrt(f);
        sJ = sqrtm(J);
        U1 = blkdiag(sJ, conj(sJ));
        V1 = blkdiag(conj(sJ), sJ);
        UU = U*U1;
        VV = V1*V;
        UU = UU(:,1:d/2);
        VV = VV(1:d/2,:);
        inj1 = real(UU);
        inj2 = imag(UU);
        prj1 = real(VV);
        prj2 = imag(VV);
        sub1 = sub.parent.subRep(inj1, 'projection', prj1, 'isIrreducible', true, 'frobeniusSchurIndicator', 1, 'isDivisionAlgebraCanonical', true);
        sub2 = sub.parent.subRep(inj2, 'projection', prj2, 'isIrreducible', true, 'frobeniusSchurIndicator', 1, 'isDivisionAlgebraCanonical', true);
        subs = {sub1 sub1};
    else % lambda < 0
         % quaternion-type representation
        assert(mod(d, 4) == 0);
        [U1, D1] = eig(J, 'vector');
        [~, I1] = sort(real(D1));
        I1 = I1(:)';
        % assert that we have (lambda, -lambda) pairs
        assert(all(D1(I1) + D1(fliplr(I1)) < tol));
        % order the pairs
        I1 = [I1(1:d/4) fliplr(I1(d/4+1:end))];
        % reorder everything and compute the left eigenvectors
        D1 = D1(I1);
        U1 = U1(:,I1);
        V1 = inv(U1);
        % chop it off
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
        W = diag([w 1i*w conj(w) -1i*conj(w)]);
        s1 = [1 0; 0 1];
        sy = [0 -1i; 1i 0];
        p = [1:d/4
             d/4+1:d/2
             d/2+1:3*d/4
             3*d/4+1:d];
        P = sparse(p(:)', 1:d, ones(1, d), d, d);
        W1 = kron(eye(d/4), [s1-sy s1+sy; -s1-sy s1-sy]/2);
        inj = U*UU * W*P*W1';
        prj = W1*P'*conj(W) * VV*V;
        assert(norm(imag(inj)) <= tol);
        assert(norm(imag(prj)) <= tol);
        inj = real(inj)/sqrt(2);
        prj = real(prj)/sqrt(2);
        subs = {sub.parent.subRep(inj, 'projection', prj, 'isIrreducible', true, 'frobeniusSchurIndicator', -2, 'isDivisionAlgebraCanonical', true)};
    end
end
