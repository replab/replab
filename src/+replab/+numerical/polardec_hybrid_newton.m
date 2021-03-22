function [U, H] = polardec_hybrid_newton(A, delta, lambda, theta)
% Polar decomposition using a Newton-type iteration
%
% Implements the polar decomposition according to
% N. J. Higham and R. S. Schreiber, "Fast Polar Decomposition of an Arbitrary Matrix"
% SIAM J. Sci. and Stat. Comput., vol. 11, no. 4, pp. 648–655, Jul. 1990
% doi: 10.1137/0911038.
%
% This is always slower than the SVD method

    d = size(A, 1);
    assert(size(A, 2) == d, 'Matrix must be square');
    if nargin < 4
        theta = 0.6;
    end
    if nargin < 3
        lambda = 0.75;
    end
    if nargin < 2
        delta = sqrt(d) * eps;
    end
    assert(delta < 1, 'Delta must be much smaller than 1.');
    X = A;
    k = 0;
    mu = 1;
    switched = false
    while mu > delta && k < 20
        k = k + 1;
        if switched
            R = eye(d) - X'*X;
            mu = norm(R, 1);
            % evaluate (3.3)
            X = X + X*R/2;
        else
            mu = normest1(@fun);
            fprintf('condest(R) = %f\n', mu);
            if mu > lambda * theta
                % evaluate (3.2)
                Xinv = inv(X);
                num = norm(Xinv,1)*norm(Xinv,Inf);
                den = norm(X,1)*norm(X,Inf);
                gamma = ((norm(Xinv, 1)*norm(Xinv, Inf))/(norm(X, 1)*norm(X, Inf)))^(1/4);
                X = (gamma*X + Xinv'/gamma)/2;
            else
                R = eye(d) - X'*X;
                mu = norm(R, 1);
                if mu > theta
                    Xinv = inv(X);
                    gamma = ((norm(Xinv, 1)*norm(Xinv, Inf))/(norm(X, 1)*norm(X, Inf)))^(1/4);
                    X = (gamma*X + Xinv'/gamma)/2;
                else
                    % evaluate (3.3)
                    X = X + X*R/2;
                    switched = true;
                end
            end
        end
    end
    U = X;
    H = (U'*A + A'*U)/2;
    function res = fun(flag, v)
        switch flag
          case 'dim'
            res = d;
          case 'real'
            res = isreal(X);
          case 'notransp'
            res = v - X'*(X*v);
          case 'transp'
            % (I - X'*X).'
            % (I - X.'*conj(X))
            res = v - X.'*(conj(X)*v);
        end
    end
end
