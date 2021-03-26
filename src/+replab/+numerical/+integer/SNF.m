function [Pcell, Qcell, IVF] = SNF(A)
%
%  This program determines the Smith Normal Form S of the M x N matrix A,
%  where S = P*A*Q with P and Q unimodular.  The invarient factors (i.e.
%  the diagonal entries of S) are returned in the array IVF.  Rather than
%  computing the matrices P and Q directly, a multipler system for P and Q
%  are returned.  This is done to avoid possible overflow (in the sense of
%  integers not being stored exactly as double precision variables) in the
%  entries of P and Q.
%
%  The multiplier system is stored in the cell arrays Pcell and Qcell.
%  The matrices P and Q can be recovered using the following:
%
%    L = min(M,N);
%    P = eye(M);
%    Q = eye(N);
%    for j=1:L
%      P = Pcell{j} * P;
%      Q = Q * Qcell{j};
%    end
%
%  Original Fortran version written by Morris Newman.
%  Adapted to Matlab with multiplier systems by Greg Wilson.
%
% See https://www.researchgate.net/post/Is-there-any-computer-algorithm-that-finds-Smith-Normal-Form-of-a-given-matrix
%
%
[M,N] = size(A);
MPLUSN=M+N;
B = zeros(MPLUSN);
B(1:M,1:N) = A;

L = min(M,N);
Pcell = cell(L,1);
Qcell = cell(L,1);
IVF = zeros(L,1);

% Set P = I and Q = I
B(1:M,N+1:MPLUSN) = eye(M);
B(M+1:MPLUSN,1:N) = eye(N);
for K=1:L
  Pcell{K} = eye(L);
  Qcell{K} = eye(L);
end

K=1;
while (K <= L)
  [I,J,T] = NMZ(K,M,N,B);
  if (T == 0)
    break;
  end

% Interchange rows i and k, if i and k are different
  if (I ~= K)
    for V=K:MPLUSN
      S      = B(K,V);
      B(K,V) = B(I,V);
      B(I,V) = S;
    end
  end

% Interchange columns j and k, if j and k are different
  if (J ~= K)
    for U=K:MPLUSN
      S      = B(U,K);
      B(U,K) = B(U,J);
      B(U,J) = S;
    end
  end

% Replace row k by its negative, if B(k,k) < 0
  if (B(K,K) < 0)
    for V=K:MPLUSN
      B(K,V) = -B(K,V);
    end
  end

% Replace row i by row i - x row k,
% x = [B(i,k)/B(k,k)], i = k+1, ... ,m

  for I=K:M
    if (I ~= K)
      X = floor(B(I,K)/T);
      for J=K:MPLUSN
        B(I,J)=B(I,J)-X*B(K,J);
      end
    end
  end

% Replace col j by col j - y col k,
% where y =[B(k,j)/B(k,k)], j = k+1, ... ,n

  U = 1;
  while (U ~= 0)
    for J=K+1:N
      Y=floor(B(K,J)/T);
      for I=K:MPLUSN
        B(I,J)=B(I,J)-Y*B(I,K);
      end
    end

%   Check to see whether or not all B(i,k) = 0,
%   i = k+1, ... ,m

    GOTO5 = 0;
    for I=K+1:M
      if (B(I,K) ~= 0); GOTO5 = 1; end
    end
    if (GOTO5 == 1); break; end

%   Check to see whether or not all B(k,j) = 0,
%   j = k+1, ... ,n

    for J=K+1:N
      if (B(K,J) ~= 0); GOTO5 = 1; end
    end
    if (GOTO5 == 1); break; end
%
%   Replace row k by row k + row p  if U <> 0
    [U] = TFD(K,M,N,B);
    if (U == 0)
      break;
    else
      for J=K:MPLUSN
        B(K,J)=B(K,J)+B(U,J);
      end
    end

  end  % end while (U ~= 0)

  if (GOTO5 == 0)
    Pcell{K} = round(B(1:M,N+1:MPLUSN));
    Qcell{K} = round(B(M+1:MPLUSN,1:N));
    B(1:M,N+1:MPLUSN) = eye(M);
    B(M+1:MPLUSN,1:N) = eye(N);
    K = K + 1;
  end

end  % end while (K <= L)

B = round(B);

% Store invariant factors in IVF
for I=1:L
  IVF(I)=B(I,I);
end

% Retrieve left and right unimodular multipliers and check result
% P = eye(M);
% Q = eye(N);
% for j=1:L
%  P = Pcell{j} * P;
%  Q = Q * Qcell{j};
% end
% S = zeros(M,N);
% for K = 1:L
%   S(K,K) = IVF(K);
% end
% norm(S - P*A*Q)

end


function [I,J,T] = NMZ(K,M,N,A)
%
%  T is set equal to 0 if Ak = 0; otherwise t is set equal to abs(A(i,j)),
%  where A(i,j) is a nonzero element which is least in absolute value.
%
for I=K:M
  for J=K:N
    T = abs(A(I,J));
    if (T ~= 0); break; end
  end
  if (T ~= 0); break; end
end
if (T == 0); return; end
for U=K:M
  for V=K:N
    S = abs(A(U,V));
    if (S ~= 0) && (S < T)
      T = S;
      I = U;
      J = V;
    end
  end
end
end


function [T] = TFD(K,M,N,A)
%
%  T is set equal to 0 if Ak is congruent to 0 modulo A(k,k); otherwise,
%  Tis set equal to p, where A(p,q) is an element of Ak which is not
%  divisible by A(k,k)
%
T=0;
for I=K:M
  for J=K:N
    if (mod(A(I,J),A(K,K)) ~= 0)
      T = I;
      return
    end
  end
end
end
