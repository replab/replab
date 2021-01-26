% Estimates the error on the coefficients of an orthogonal matrix sampled from the Haar measure
% using the technique in http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
%
% We sample a real matrix from the Haar measure, orthogonal up to double floating point precision,
% and then use https://en.wikipedia.org/wiki/Orthogonal_matrix#Nearest_orthogonal_matrix to compute
% an orthgonal matrix in GEM default precision (50 digits)

addpath('~/software/gem/gem');
values_n = [2 4 8 16 32 64]; % use powers of two, so that nSamples is integer
errors = [];
nCoeffs = values_n(end)^2;
for i = 1:length(values_n)
    n = values_n(i);
    nSamples = nCoeffs / (n*n); % we want the same number of errors for each size
    err = [];
    frob_err = [];
    fprintf('Computing for n=%d\n', n);
    for j = 1:nSamples
        % sample
        X = randn(n, n);
        [Q, R] = qr(X);
        F = diag(diag(R)./abs(diag(R)));
        Q = Q*F;
        % save original sample in double precision
        Q0 = Q;
        % Newton orthogonalization iterations
        Q = gem(Q);
        nIters = 3;
        for k = 1:nIters
            N = Q'*Q;
            P = Q*N/2;
            Q = 2*Q + P*N - 3*P;
        end
        errQ = double(abs(Q-Q0)); % compute the error
        err = [err; errQ(:)]; % tabulate
    end
    errors = [errors err(:)];
end
save -v7 orthogonalError.mat errors values_n