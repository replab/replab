function [d, X] = xgcd(A)
% Returns the greatest common divisor of an integer array and the coefficients of Bezout's identity
%
% A better approach with smaller coefficients would be given in:
% B. S. Majewski and G. Havas, "The complexity of greatest common divisor computations"
% in Algorithmic Number Theory, Berlin, Heidelberg, 1994, pp. 184–193, doi: 10.1007/3-540-58691-1_56.
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
        [d, x1, x2] = gcd(A(1), A(2));
        X = [x1 x2];
      otherwise
        [dr, Xr] = replab.numerical.integer.xgcd(A(2:end));
        [d, x1, x2] = gcd(A(1), dr);
        X = [x1 x2*Xr];
    end
end
