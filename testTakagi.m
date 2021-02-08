n = 200;
X = randn(n,n) + 1i*randn(n,n);
J = X.'+X;
tic;
[U,D] = replab.numerical.takagi(J, 'algorithm', 'svd');
toc
norm(U.'*D*U-J, 'fro')
tic;
[U,D] = replab.numerical.takagi(J, 'algorithm', 'jacobi');
toc
norm(U.'*D*U-J, 'fro')