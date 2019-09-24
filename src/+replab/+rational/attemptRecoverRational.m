function [num den] = attemptRecoverRational(M)
% Attempts to recover a rational approximation of a double vector/matrix
%
% Args:
%  M (double vector or matrix): Vector or matrix to find a rational approximation of
%
% Returns
% -------
%   num: double vector or matrix
%     Numerator of the approximation with the same shape as `M`
%   den: double scalar
%     Denominator such that ``num/den`` approximates `M`

    tol = replab.Settings.doubleEigTol; % magic epsilon to remove
    maxden = 1200; % maximum denominator

    v = M(:); % vectorize the input
    d = length(v);
    numv = zeros(d, 1); % numerators
    denv = zeros(d, 1); % denominators
    [v1, I] = sort(v); % sort the coefficients to group them when they are close
    start = 1;
    untl = 1;
    den = 1; % current common denominator
    while start <= d
        % finds a run of values that are closed to each other
        while untl <= d && v1(untl) - v1(start) < tol
            untl = untl + 1;
        end
        range = start:(untl-1);
        m = mean(v1(range)); % take the average to reduce the noise
        [num1 den1] = rat(m, tol); % find a rational approximation
        den = lcm(den, den1); % update the common denominator
        if den > maxden
            % common denominator too big? bail out
            num = [];
            den = [];
            return
        end
        % update rational decomposition in the original ordering
        numv(I(range)) = num1;
        denv(I(range)) = den1;
        % start the process with the new block
        start = untl;
    end
    num = reshape(numv .* (den./denv), size(M));
end
