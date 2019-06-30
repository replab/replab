function [V, D] = eig(A)
% EIG    Eigenvalues and eigenvectors.
% (Quaternion overloading of standard Matlab function, with limitations.)
%
% Acceptable calling sequences are: [V,D] = EIG(X) and V = EIG(X).
% The results are as for the standard Matlab EIG function (q.v.).

% Copyright (c) 2005 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 2)

[r, c] = size(A);

if r ~= c
    error('Matrix must be square');
end

if r == 1
    
    % The argument is a single quaternion. This case could be handled by
    % using the standard Matlab eig() function on the complex adjoint of
    % A, but there are problems if A is a complexified quaternion, since
    % we cannot make a complex value with complex parts.
    %
    % For this reason we output an error message and leave it to the user
    % to use the appropriate adjoint.
    
    A = inputname(1); if A == ''; A = 'A'; end
    
    disp('The eig function is not implemented for single quaternions.');
    disp(sprintf('Try using eig(adjoint(%s, ''real'')) or', A));
    disp(sprintf('          eig(adjoint(%s, ''complex'')',  A));
    error('Implementation restriction - see advice above');
    
    % TODO The paper below, Theorem 4.7, may offer a way to compute a
    % decomposition of a single quaternion. To be studied. Whether it
    % belongs in this function is to be considered, it could be a separate
    % function. A key issue is whether Oba's decomposition would match with
    % one computed using a matrix representation of the quaternion,
    % e.g computed by the QTFM function ADJOINT.
    
    % Roger M. Oba, 'Eigen-Decomposition of Quaternions',
    % Advances in Applied Clifford Algebras, Oct 2018, 28(5), p94.
    % doi:10.1007/s00006-018-0911-6
end

if ~isreal(A)
   warning(['EIG does not give correct results for a complex ',...
            'quaternion  matrix. See TODO note in code.'])
   % TODO
   %
   % The problem is with the definition of 'Hermitian', and the eigenvalue
   % decomposition. The following article may be of relevance to this issue,
   % and it may require the functon tridiagonalize to be modified.
   %
   % Yongge Tian, Matrix Theory over the Complex Quaternion Algebra,
   % arXiv:math/0004005v1 [math.RA], 1 Apr 2000.
end

% The method used depends on whether A is Hermitian or otherwise.

if ishermitian(A)
    
    % For a Hermitian matrix the eigenvalues and eigenvectors are real if A is a
    % real quaternion matrix and complex if A is a complexified quaternion matrix.
    % The first step is to tridiagonalize A using Householder transformations to
    % obtain P and B where B is a real or complex tridiagonal symmetric matrix,
    % and P is the product of the Householder matrices used to compute B
    % such that P'*B*P = A.

    [P, B] = tridiagonalize(A);
else
    
    % For a general (i.e. non Hermitian) matrix, we don't currently have a method to
    % convert A to a real matrix. The only feasible approach is to use an adjoint
    % matrix, but we leave this to the user.
    
    disp('The eig function is not implemented for non-Hermitian quaternion matrices.');
    disp('Eigenvalues and eigenvectors can be computed by using the standard Matlab');
    disp('eig function on an adjoint matrix by using the function adjoint() (q.v.).');
    error('Implementation restriction - see advice above');

end

% We now have: P' * B * P = A.

% The second step is to compute the eigenvectors and eigenvalues of B using the standard
% Matlab routine, for a real or complex matrix (which happens to be tridiagonal).

if nargout == 0
    eig(B)
elseif nargout == 1
    V = eig(B);
else    
    [V, D] = eig(B);
    V = P' * V; % Combine P (from the tridiagonalization) with V.
end

% TODO
%
% Consider:
%
% F. O. Farid , Qing-Wen Wang & Fuzhen Zhang (2011): On the eigenvalues of
% quaternion matrices, Linear and Multilinear Algebra, 59:4, 451-473
% DOI:10.1080/03081081003739204

% $Id: eig.m 1016 2019-01-17 10:35:20Z sangwine $
