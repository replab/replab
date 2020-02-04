function err = errorModel(X)
% Rough error model used for computations that have no approximation error
%
% This is supposed to represent floating point approximation error. In practice,
% only (signed) permutation representations are practical in RepLAB, so this
% overestimates the error.
    geod = sqrt(size(X, 1) * size(X, 2));
    err = norm(X, 'fro') * dav * 5e-16;
end
