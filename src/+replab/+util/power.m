function res = power(a, n, op)
% Computes ``a^n`` by repeated squaring
%
% Args:
%   a: Element to compute the power of
%   n (integer): Power, must be > 0
%   op (function_handle): Binary operation representing multiplication
    if n == 1
        res = a;
    else
        b = a;
        extra = a;
        k = n - 1;
        while k > 1
            if mod(k, 2) == 1
                x = op(b, extra);
                k = (k-1)/2;
            else
                x = extra;
                k = k/2;
            end
            b = op(b, b);
            extra = x;
        end
        res = op(b, extra);
    end
end
