function [sub1, sub2] = regularizeRealPair(sub)
% Splits a subrepresentation encoding, using a complex division algebra, a pair of real subrepresentations
%
% Args:
%   sub (`+replab.SubRep`): Subrepresentation with its ``divisionAlgebraName`` set to ``'complex'``
%
% Returns
% -------
%   sub1: `+replab.SubRep`
%     First real-type irreducible subrepresentation
%   sub2: `+replab.SubRep`
%     Second real-type irreducible subrepresentation
    X = sub.parent.commutant.sample;
    d = sub.dimension;
    I = sub.injection;
    P = sub.projection;
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
    s1 = sqrt(K1);
    s3 = sqrt(K3);
    Tev = diag([conj(s1);inv(s3);conj(s3)]);
    Tiev = diag([s1;s3;inv(conj(s3))]);
    % put in diagonal form
    C = [1+1i 1-1i; 1-1i 1+1i]/2;
    Ci = conj(C);
    Tform = blkdiag(eye(n1), kron(C, eye(n23)));
    Tiform = blkdiag(eye(n1), kron(Ci, eye(n23)));
    T = blkdiag(conj(Tev*Tform), Tev*Tform);
    Ti = blkdiag(conj(Tiform*Tiev), Tiform*Tiev);
    inj = U*W*T*kron([1 1i; 1 -1i], eye(d/2))/sqrt(2);
    prj = kron([1 1; -1i 1i], eye(d/2))*Ti*Wi*Ui/sqrt(2);
    tol = 1e-10;
    assert(norm(imag(inj)) < tol);
    assert(norm(imag(prj)) < tol);
    inj = real(inj);
    prj = real(prj);
    sub1 = sub.parent.subRep(inj(:,1:d/2), 'projection', prj(1:d/2,:), 'isIrreducible', true, 'frobeniusSchurIndicator', 1);
    sub2 = sub.parent.subRep(inj(:,d/2+1:end), 'projection', prj(d/2+1:end,:), 'isIrreducible', true, 'frobeniusSchurIndicator', 1);
end
