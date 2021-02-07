function s1 = regularizeQuaternion(s)
% Regularizes the given quaternion representation, whose complex algebra has been revealed
    X = s.parent.commutant.sample;
    X1 = s.parent.commutant.sample;
    d = s.dimension;
    I = s.injection;
    P = s.projection;
    % get the real and imaginary parts for both the injection and projection maps
    A = I(:,1:2:d);
    B = I(:,2:2:d);
    C = P(1:2:d,:);
    D = P(2:2:d,:);
    % construct the complex change of basis
    U = [A+1i*B A-1i*B]/sqrt(2);
    Ui = [C-1i*D; C+1i*D]/sqrt(2);
    % and express the sample in that basis
    Y = Ui*X*U;
    % the sample has the shape, in that basis
    % [x*I     J
    %  conj(J) conj(x)*I] where I is the identity matrix of the relevant size
    J = Y(1:d/2,d/2+1:end);
    % and because this is a quaternionic representation, we have conj(J)*J = -f with f > 0
    f = -trace(conj(J)*J)/(d/2);
    % so we normalize J
    J = J / sqrt(f);
    % thus we have conj(J)*J = -I
    %
    % now, let v_i be an eigenvector of J with eigenvalue l_i: J*v_i = l_i*v_i
    %
    % then we have conj(J)*v_i = -l_i*v_i
    % which is also written J*conj(v_i) = -conj(l_i)*v_i
    % and thus conj(v_i) is also an eigenvector of J with eigenvalue -1/conj(l_i)
    %
    % we do the EV decomposition of J
    [V,ev] = eig(J);
    Vi = inv(V);
    % we remark that Y has the shape
    % [x*I          V*ev*Vi
    %  conj(V*ev*Vi) conj(x)*I]
    ev = diag(ev);
    [~, I1] = sort(real(ev));
    [~, I2] = sort(-real(1./conj(ev)));
    assert(all(I1 ~= I2));
    mask = I1 < I2;
    I = [I1(mask);I2(mask)];
    V = V(:,I);
    Vi = inv(V);
    W = blkdiag(conj(V),V);
    Wi = blkdiag(conj(Vi),Vi);
    Y1 =  Wi*Ui*X*U*W/sqrt(f);
    L = Y1(1:d/4,d/2+d/4+1:end);
    K = Y1(d/4+1:d/2,d/2+1:d/2+d/4);
    L = diag(diag(L));
    T = blkdiag(L, eye(d/4), conj(L), eye(d/4));
    Ti = blkdiag(inv(L), eye(d/4), conj(inv(L)), eye(d/4));
    % now Wi*Ui*X*U*W/sqrt(f) has the shape
    %
    % [    x*I       0         0         L
    %        0     x*I         K         0
    %        0 conj(L) conj(x)*I         0
    %  conj(K)       0         0 conj(x)*I  ]
    %
    % and we are making that matrix look like
    %
    % [      x*I          0         0       y*I
    %          0        x*I      -y*I         0
    %          0 conj(-y*I) conj(x)*I         0
    %  conj(y*I)          0         0 conj(x)*I  ]
    %
    %
    % computing inv(diag(E,F,conj(E),conj(F))) * Z * diag(E,F,conj(E),conj(F))
    % with E,F diagonal matrices, we get
    %
    % [                  inv(E)*L*conj(F)
    %              inv(F)*K*conj(E)
    %         conj(inv(E)*L)*F
    %   conj(inv(F)*K)*E                  ]  + blkdiag(x*I, x*I, conj(x)*I, conj(x)*I)
    %
    % so we set F = I, E = L
    s1 = [1 0; 0 1];
    sy = [0 -1i; 1i 0];
    S = 1i*kron([s1-sy 1i*(s1+sy); -s1-sy 1i*(s1-sy)], eye(d/4))/2;
    S = reshape(permute(reshape(S, [d 4 d/4]), [1 3 2]), [d d]);
    Si = S';
    Y = Si*Ti*Wi*Ui*X*U*W*T*S;
    inj = U*W*T*S;
    prj = Si*Ti*Wi*Ui;
    tol = 1e-10;
    assert(norm(imag(inj)) < tol);
    assert(norm(imag(prj)) < tol);
    inj = real(inj);
    prj = real(prj);
    s1 = s.parent.subRep(inj, 'projection', prj, 'divisionAlgebraName', 'quaternion.rep', 'isIrreducible', true, 'frobeniusSchurIndicator', -2);
end
