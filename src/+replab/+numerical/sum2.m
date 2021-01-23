function [res, err] = sum2(p)
% Compensated summation of a vector
%
% Implements the algorithm Sum2 from
% T. Ogita, S.M. Rump and S. Oishi, Accurate Sum and Dot Product
% SIAM Journal on Scientific Computing, 26(6):1955-1988, 2005
% https://doi.org/10.1137/030601818
%
% Args:
%   p (double(1,\*)): Real or complex vector of coefficients to sum
%
% Returns
% -------
%   res: double
%     Result of summation
%   err: double
%     Residual error

    e = 1e-20;
    if 1+e ~= 1-e
        warning('Rounding mode is not round to nearest, results may be inaccurate');
    end

    acc1 = p(1);
    acc2 = 0;
    for i = 2:length(p)
        x = acc1 + p(i); % TwoSum(acc1, p(i))
        z = x - acc1;
        y = (acc1 - (x - z)) + (p(i) - z);
        % result of TwoSum in [x, y]
        acc2 = acc2 + y;
        acc1 = x;
    end
    % instead of returning just acc1+acc2, we perform TwoSum to return the error as well
    x = acc1 + acc2; % TwoSum(acc1, acc2)
    z = x - acc1;
    y = (acc1 - (x - z)) + (acc2 - z);
    res = x; % result
    err = y;
end
