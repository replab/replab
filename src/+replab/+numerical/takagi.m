function [U, D] = takagi(A, varargin)
% Returns the Takagi factorization of the given complex symmetric matrix
%
% We propose two methods:
%
% Keyword Args:
%   algorithm ('svd', 'jacobi', optional): Algorithm, default: svd
%   maxSweeps (integer, optional): Jacobi algorithm, maximal number of sweeps across the matrix, default: 50
%
% A = U.'*D*U
    n = size(A, 1);
    assert(size(A, 2) == n, 'A must be square');
    if n == 1
        U = conj(sqrt(sign(A)));
        D = abs(A);
        return
    end
    assert(n >= 2, 'A must be at least 2x2');
    assert(all(all(A == A.'))); % A is symmetric
    args = struct('algorithm', 'svd', 'maxSweeps', 50);
    args = replab.util.populateStruct(args, varargin);
    maxSweeps = args.maxSweeps;
    algorithm = args.algorithm;
    switch algorithm
      case 'svd'
        [U,S,V] = svd(A);
        Z = U'*conj(V);
        Uz = U*sqrtm(Z);
        U = Uz.';
        d = diag(S);
      case 'jacobi'

        d = zeros(1, n);
        ev = zeros(2, n);
        ev(2, :) = diag(A);
        U = eye(n);

        red = 0.04/n^4;

        done = false;

        sym_eps = 2*eps^2;
        sq_eps = 4*eps^2;

        nsweeps = 1;

        while ~done && nsweeps <= maxSweeps

            T = triu(A, 1); % upper triangle
            off = sum(sum(conj(T).*T)); % sum(abs2(...)
            if off <= sym_eps
                done = true;
                continue
            end

            if nsweeps < 4
                thresh = off*red;
            else
                thresh = 0;
            end

            for q = 2:n
                for p = 1:q-1
                    off = conj(A(p,q))*A(p,q);
                    sqp = conj(ev(2,p))*ev(2,p);
                    sqq = conj(ev(2,q))*ev(2,q);
                    if nsweeps > 4 && off < sym_eps*(sqp+sqq)
                        A(p,q) = 0;
                    elseif off > thresh
                        t = abs(sqp-sqq) / 2;
                        if t > eps
                            f = sign(sqp-sqq) * (ev(2,q)*conj(A(p,q)) + conj(ev(2,p))*A(p,q));
                        else
                            if sqp == 0
                                f = 1;
                            else
                                f = sqrt(ev(2,q)/ev(2,p));
                            end
                        end

                        t = t + sqrt(t^2 + conj(f)*f);
                        f = f / t;

                        ev(1,p) = ev(1,p) + A(p,q)*conj(f);
                        ev(2,p) = A(p,p) + ev(1,p);
                        ev(1,q) = ev(1,q) - A(p,q)*f;
                        ev(2,q) = A(q,q) + ev(1,q);

                        t = conj(f)*f;
                        ci = sqrt(t + 1);
                        f = f / ci;
                        t = t / (ci*(ci+1));

                        x = A(1:p-1,p);
                        y = A(1:p-1,q);
                        A(1:p-1,p) = x + (conj(f)*y - t*x);
                        A(1:p-1,q) = y - (f*x + t*y);

                        x = A(p,p+1:q-1);
                        y = A(p+1:q-1,q).';
                        A(p,p+1:q-1) = x + (conj(f)*y - t*x);
                        A(p+1:q-1,q) = y - (f*x + t*y);

                        x = A(p,q+1:n);
                        y = A(q,q+1:n);
                        A(p,q+1:n) = x + (conj(f)*y - t*x);
                        A(q,q+1:n) = y - (f*x + t*y);

                        A(p,q) = 0;

                        x = U(p,:);
                        y = U(q,:);
                        U(p,:) = x + (f*y - t*x);
                        U(q,:) = y - (conj(f)*x + t*y);
                    end % elseif off > thresh
                end % for p=1:q
            end % for q=2:n

            for p = 1:n
                ev(1,p) = 0;
                A(p,p) = ev(2,p);
            end

            nsweeps = nsweeps + 1;
        end % while

        if ~done
            error('Bad convergence in takagi');
        end

        % Make the diagonal elements non-negative
        for p = 1:n
            d(p) = abs(A(p,p));
            if d(p) > eps && d(p) ~= real(A(p,p))
                f = sqrt(A(p,p)/d(p));
                for q = 1:n
                    U(p,q) = U(p,q) * f;
                end
            end
        end

    end
    D = diag(d);
end
