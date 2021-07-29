function [H, U] = hermiteNormalForm(A)
% Computation of the (row-style) Hermite normal form using the Cohen algorithm
%
% It returns a unimodular matrix ``U`` and the Hermite Normal Form ``H``, such that ``H == U*A``.
%
% This code is taken from AbstractAlgebra.jl `<https://github.com/wbhart/AbstractAlgebra.jl/>`_ under
% the Simplified "2-clause" BSD License (as its individual files are).
%
% The code has been translated from the original Julia code.
%
% Example:
%   >>> [H, U] = replab.numerical.integer.hermiteNormalForm([3 3 1 4; 0 1 0 0; 0 0 19 16; 0 0 0 3]);
%   >>> all(all(H == [3 0 1 1; 0 1 0 0; 0 0 19 1; 0 0 0 3]))
%       1
%   >>> [H, U] = replab.numerical.integer.hermiteNormalForm([2 3 6 2; 5 6 1 6; 8 3 1 1]);
%   >>> all(all(H == [1 0 50 -11; 0 3 28 -2; 0 0 61 -13]))
%       1
%
% Args:
%   integer(m,n): Integer matrix to decompose
%
% Returns
% -------
%   H: integer(m,n)
%     Hermite Normal Form
%   U: integer(m,m)
%     Unimodular matrix

    if round(A) ~= A
        error('Requires integer input.');
    end

    H = A;
    m = size(H, 1);
    n = size(H, 2);
    U = eye(m); % The Julia code admits U as a parameter, otherwise initialized to the identity
    l = min(m, n);
    k = 1;
    t = [];
    t1 = [];
    t2 = [];
    for i = 1:l
        for j = k + 1:m
            if H(j, i) == 0
                continue
            end
            [d, u, v] = gcd(H(k, i), H(j, i));
            a = H(k, i)/d;
            b = -H(j, i)/d;
            for c = i:n
                t = H(j, c);
                t1 = a * H(j, c);
                t2 = b * H(k, c);
                H(j, c) = t1 + t2;
                t1 = u * H(k, c);
                t2 = v * t;
                H(k, c) = t1 + t2;
            end
            for c = 1:m
                t = U(j, c);
                t1 = a * U(j, c);
                t2 = b * U(k, c);
                U(j, c) = t1 + t2;
                t1 = u * U(k, c);
                t2 = v * t;
                U(k, c) = t1 + t2;
            end
        end
        if H(k, i) == 0
            continue
        end
        cu = sign(H(k, i)); % canonical_unit
        if cu ~= 1
            for c = i:n
                H(k, c) = H(k, c) / cu;
            end
            for c = 1:m
                U(k, c) = U(k, c) / cu;
            end
        end
        for j = 1:k-1
            q = -floor(H(j, i)/H(k, i));
            for c = i:n
                t = q * H(k, c);
                H(j, c) = H(j, c) + t;
            end
            for c = 1:m
                t = q * U(k, c);
                U(j, c) = U(j, c) + t;
            end
        end
        k = k + 1;
    end
end