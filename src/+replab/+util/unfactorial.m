function [n r] = unfactorial(f)
% Computes the inverse factorial
%
% Returns $n$ and $r$ such that $n$ is the largest integer with $f = n! + r$.
%
% Args:
%   f (vpi): Number to compute the factorial of
%
% Returns
% -------
%   n:
%     integer: Factorial preimage
%   r:
%     vpi or double: Remainder (of same type as ``f``)
    n = 1;
    m = 1;
    if isa(f, 'vpi');
        n = vpi(n);
        m = vpi(m);
    end
    while 1
        newM = (n + 1) * m;
        if newM > f
            r = f - m;
            return
        end
        m = newM;
        n = n + 1;
    end
    n = double(n);
end
