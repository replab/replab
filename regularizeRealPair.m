function s1 = regularizeRealPair(s)
% Regularizes the given pair of real representations, whose complex algebra has been revealed
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
    % and because this is a pair of real representations, we have conj(J)*J = f with f > 0
    f = trace(conj(J)*J)/(d/2);
    % so we normalize J
    J = J / sqrt(f);
    % thus we have conj(J)*J = I
    %
    % now, let v_i be an eigenvector of J with eigenvalue l_i: J*v_i = l_i*v_i
    % then we have conj(J)*v_i = l_i*v_i
    % which is also written J*conj(v_i) = conj(l_i)*v_i
    % and thus conj(v_i) is also an eigenvector of J with eigenvalue 1/conj(l_i)
    %
    % there are two cases:
    % - either v_i is real and thus v_i=conj(v_i); then abs(l_i) = 1
    % - either v_i is complex, and there is another conj(v_i) eigenvector
    %
    % we do the EV decomposition of J
    [V,ev] = eig(J);
    Vi = inv(V);
    ev = diag(ev);
    [~, ind1] = sort(real(ev));
    [~, ind2] = sort(real(1./conj(ev)));
    I1 = ind1(ind1 == ind2);
    I2 = ind1(ind1 < ind2);
    I3 = ind1(ind1 > ind2);
    n1 = length(I1);
    n23 = length(I2);
    I = [I1; I2; I3];
    V = V(:,I);
    Vi = inv(V);
    W = blkdiag(conj(V),V);
    Wi = blkdiag(conj(Vi),Vi);
    Y =  Wi*Ui*X*U*W/sqrt(f);
    K = Y(1:d/2, d/2+1:end);
    % K has the form
    %
    % [K1  0  0
    %   0  0 K2
    %   0 K3  0] where K1,K2,K3 are diagonal matrices
    K1 = diag(K(1:n1,1:n1));
    K2 = diag(K(n1+1:n1+n23,n1+n23+1:end));
    K3 = diag(K(n1+n23+1:end,n1+1:n1+n23));
    % correct eigenvalues
    Tev = diag([conj(sqrt(K1));inv(sqrt(K3));conj(sqrt(K3))]);
    % put in diagonal form
    C =
    Tform = blkdiag(eye(n1), conj(kron([-1+1i -1-1i; -1-1i -1+1i], eye(n23))));
    T = blkdiag(conj(Tev*Tform), Tev*Tform);
    inv(T)*Y*T
    error('asd')
    %y12 = sqrt(diag(Y12));
    %T = diag([y12(1:length(I1)); ones(d/2-length(I1), 1); conj(y12(1:length(I1))); ones(d/2-length(I1), 1)]);
    %Y =  inv(T)*Wi*Ui*X*U*W*T/sqrt(f);
    %Y1 =  inv(T)*Wi*Ui*X1*U*W*T/sqrt(f);
end
