function [U,S,V] = svdj(A, tol)
% SVDJ Singular value decomposition using Jacobi algorithm.

% Copyright (c) 2005, 2006 Nicolas Le Bihan and Stephen J. Sangwine.

% This function works for real, complex or quaternion matrices.
%
% Arguments: A    - a real, complex or quaternion matrix.
%            tol  - a tolerance, defaults to eps if omitted.
% Returns:   U, V - singular vectors (unitary matrices).
%            S    - singular values.
%
% The singular value decomposition computed by this function is the same as
% that computed by the function svd.m, but this code, being based on a
% cyclic Jacobi algorithm, is more accurate. However, it is also slower.
% It is intended as a reference implementation since the Jacobi algorithm
% is known to be the most accurate SVD algorithm.
%
% This function will work on real, complex and quaternion matrices.
%
% References:
%
% N. Le Bihan and S. J. Sangwine, "Jacobi Method for Quaternion Matrix
% Singular Value Decomposition", Applied Mathematics and Computation,
% 187(2), 15 April 2007, 1265?1271.   DOI:10.1016/j.amc.2006.09.055.
%
% S. J. Sangwine and N. Le Bihan, "Computing the SVD of a quaternion matrix", 
% in Seventh International Conference on Mathematics in Signal Processing,
% 17-20 December 2006, The Royal Agricultural College, Cirencester, UK,
% pp5-8. Institute of Mathematics and its Applications, 2006.

narginchk(1, 2), nargoutchk(0, 3)

if nargin == 1
    tol = eps; % Default value for the tolerance.
end

if nargout == 2
   error('The number of output parameters must be 0, 1 or 3'); 
end

% We should verify here that A is a real/complex (i.e. double) or
% quaternion matrix. The code cannot work for arbitrary datatypes.

if ~isreal(A) && isa(A,'quaternion')
   error('svd_jacobi does not work with complex quaternion matrices'); 
end

[M,N] = size(A); K = min(M,N); % K is the number of singular values.

% In what follows we need to be able to construct a quaternion or real or
% complex matrix according to the type of the actual supplied for the
% parameter A. This is a tricky bit of programming, but the key to it is
% the concept of function handles (see Matlab Help). This permits us to
% call the constructor function for the type of A, which is double in the
% case of real or complex A, quaternion if A is quaternion.

F = str2func(class(A)); % F is a function handle.

V = F(eye(N));

% Calculate the sum of the moduli of the diagonal elements of the
% implicit matrix B = A' * A. This sum is invariant and we do not
% need to calculate it again. We normalise it by the matrix size.

On = 0; for c = A, On = On + sum(abs(c).^2); end; On = On ./ N;

Previous_Off = Inf; % We test on each sweep to make sure the Off diagonal
                    % sum is reducing. If it does not reduce we stop. We
                    % use infinity as the initial value.
  
while true % This is really a repeat .. until, but since Matlab does not
           % provide this statement, we use an if .. break at the end of
           % the loop.

  % Sweep through the upper triangular part of the implicit matrix B.
  
  R = 0; % We count the rotations so that we know if we have not done any
         % during a whole sweep.
  
  for r = 1 : N - 1
    for c = r + 1 : N

      % Calculate the three elements of the implicit matrix B that are
      % needed to calculate a Jacobi rotation. Since B is Hermitian, the
      % fourth element (b_cr) is not needed.
      
      b_rr = sum(abs(A(:,r)).^2); % Real value.
      b_cc = sum(abs(A(:,c)).^2); % Real value.
      b_rc = A(:,r)' * A(:,c);    % Same type as A.

      % Calculate a Jacobi rotation (four elements of G).  The two values
      % that we calculate are a real value, C = cos(theta) and S, a value
      % of the same type as A, such that |S| = sin(theta).
      
	  m = abs(b_rc);
      
      if m ~= 0 % If the off-diagonal element is zero, we don't rotate.

        tau = (b_cc - b_rr) / (2 * m); % tau is real and will be zero if
                                       % the two on-diagonal elements are
                                       % equal. In this case G will be an
                                       % identity matrix, and there is no
                                       % point in further calculating it.
        if tau ~= 0
          
          R = R + 1; % Count the rotation we are about to perform.

		  t   = sign(tau) ./ (abs(tau) + sqrt(1 + tau .^ 2));
		  C   = 1 ./ sqrt(1 + t .^ 2); 
		  S   = (b_rc .* t .* C) ./ m;
        
          % Initialize the rotation matrix, which is the same size as the
          % implicit matrix B.

          % We have to create an identity matrix here of the same type as A,
          % that is, quaternion if A is a quaternion, double if A is double.
          % To do this we use a function handle (q.v.) constructed from the
          % class type of A. This was done before the loop, since the type
          % of A is invariant.

          G = F(eye(N));

          G(r,r) = F(C);
          G(c,c) = F(C);
          G(r,c) = S;
          G(c,r) =-conj(S);
     
          % Update of A and V. This is apparently terrible as it performs a
          % full matrix multiplication and most of G is zero. An alternative
          % is to multiply only the row and column of A that are affected,
          % but this turns out to work slower by about a factor of two. See
          % the (commented) code below. Clearly Matlab exploits the near
          % identity structure of G in some way, perhaps because it uses a
          % special (sparse?) structure to store a near identity matrix.
          %
          % Temp   = A(:,r) .* G(r,r) + A(:,c) .* G(c,r);
          % A(:,c) = A(:,r) .* G(r,c) + A(:,c) .* G(c,c);
          % A(:,r) = Temp;

          A = A * G;
          V = V * G;
        end
      end
    end
  end
  
  if R == 0 % then there were no rotations on the last sweep...
      
      % This condition can occur with pathological matrices and there must
      % be some fix to enable the SVD to be calculated. For the moment at
      % least the problem is detected - prior code would enter an infinite
      % loop under this condition.
      
      error('No rotations performed during sweep.')
  end

  % Calculate the sum of the off-diagonal elements of the matrix B.

  B = A' * A;
  
  Off = sum(sum(abs(triu(B, 1))))/(N.^2); % Normalise by the matrix size!

  if (Off/On) < tol
    break; % The off-diagonal sum is small enough for us to stop.
  end
  
  if Previous_Off < Off
      warning('QTFM:information', ...
              'Terminating sweeps: off diagonal sum increased on last sweep.')
      break;
  end;
  
  Previous_Off = Off;
  
end

% Extract and sort the singular values. The vector T may be longer than the
% number of singular values (K) in cases where A is not square.
                        
[T,IX] = sort(sqrt(abs(diag(B))),'descend');

if nargout == 0 || nargout == 1 % .. only the singular values are needed.
  U = T(1:K);
end

if nargout == 3 % .. the singular vectors and singular values are needed.

    A = A(:, IX); % Map the columns of A and V into the same order as the
    V = V(:, IX); % singular values, using the sort indices in IX.
    
    % Construct the left singular vectors. These are in A but we need
    % to divide each column by the corresponding singular value. This
    % calculation is done by replicating T to make a matrix which can
    % then be divided into A element-by-element (vectorized division).

    U = A ./ repmat(T',M,1);

    S = diag(T);  % Construct a diagonal matrix of singular values from
                  % the vector T, because when there are three output
                  % parameters, S is required to be a matrix.
end
% $Id: svdj.m 1004 2017-11-15 17:14:09Z sangwine $

