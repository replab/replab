function t0 = tinv(p0, v)
% Inverse of Student's T cumulative distribution function
%
% Based on https://www.mathworks.com/matlabcentral/answers/260564-how-to-implement-tinv-manually-without-statistics-and-machine-learning-toolbox
%
% Note: we do not accept vectorized inputs.
%
% Args:
%   p0 (double): Resulting value
%   v (integer): Number of degrees of freedom
%
% Returns:
%   double: Inverse
    if p0 < 0.5
        t0 = -replab.numerical.tinv(1 - p0, v);
    else
        t0 = fzero(@(u) p0 - cdf(u, v), 5);
    end
end

function p = cdf(t, v)
    p2 = 1 - betainc(v/(v+t^2), v/2, 1/2); % 2-tailed t-distribution
    p = 1 - (1-p2)/2; % 1-tailed t-distribution
end
