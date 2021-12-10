function Y = recoverRational(X, tol, maximumDenominator)
% Recovers a rational approximation of a floating-point matrix
%
% To save on the number of computations, this function groups clusters coefficients that are closer than the given tolerance.
%
% Args:
%   X (double(\*,\*)): Matrix to find an approximation of
%   tol (double): Tolerance
%   maximumDenominator (integer): Maximum denominator
%
% Returns:
%   `+replab.cyclotomic` (\*,\*) or ``[]``: Rational approximation if one is found
    assert(maximumDenominator < 2^53, 'Maximal denominator needs to be exactly representable in double floating-point arithmetic');
    v = X(:); % vectorize the input
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
        if max(v1(range)) - min(v1(range)) > tol
            Y = [];
            return
        end
        m = mean(v1(range)); % take the average to reduce the noise
        [num1, den1] = rat(m, tol); % find a rational approximation
        if den1 > maximumDenominator
            Y = [];
            return
        end
        % update rational decomposition in the original ordering
        numv(I(range)) = num1;
        denv(I(range)) = den1;
        % start the process with the new block
        start = untl;
    end
    numv = reshape(numv, size(X));
    denv = reshape(denv, size(X));
    Y = replab.cyclotomic.fromRationals(numv, denv);
end
