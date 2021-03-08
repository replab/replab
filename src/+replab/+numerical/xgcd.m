function [d, X] = xgcd(A)
% Returns the greatest common divisor of an integer array and the coefficients of Bezout's identity
%
% Args:
%   A (integer(1,\*)): Integer array, elements can be negative or zero as well
%
% Returns
% -------
%   d: integer
%     Greatest common divisor
%   X: integer(1,\*)
%     Coefficients such that ``dot(X, A) = d``
    switch length(A)
      case 1
        d = abs(A);
        X = sign(A);
      case 2
        [d, x1, x2] = xgcd2(A(1), A(2));
        X = [x1 x2];
      otherwise
        [dr, Xr] = replab.numerical.xgcd(A(2:end));
        [d, x1, x2] = xgcd2(A(1), dr);
        X = [x1 x2*Xr];
    end
end

function [d,x,y] = xgcd2(a, b)
% Returns the greatest common divisor of two integers and the coefficients of Bezout's identity
    if b == 0
        x = 1;
        y = 0;
        d = a;
        return
    end
    [d, x1, y1]= xgcd2(b, mod(a,b));
    x = y1;
    y = x1 - floor(a/b)*y1;
    if d < 0
        x1 = -x1;
        y1 = -y1;
        d = -d;
    end
end
