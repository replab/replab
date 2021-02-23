n = 100; % matrix size

X = randn(n,n) + 1i*randn(n,n);
J = X.'+X;

disp(' ');
disp('SVD Algorithm');
disp('=============');
tic;
[U,D] = replab.numerical.takagi(J, 'algorithm', 'svd');
toc
norm(U.'*D*U-J, 'fro')

disp(' ');
disp('Hybrid algorithm');
disp('================');
tic;
[U,D] = replab.numerical.takagi(J, 'algorithm', 'hybrid');
toc
norm(U.'*D*U-J, 'fro')

disp(' ');
disp('Jacobi algorithm');
disp('================');
tic;
[U,D] = replab.numerical.takagi(J, 'algorithm', 'jacobi');
toc
norm(U.'*D*U-J, 'fro')
