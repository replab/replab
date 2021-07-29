% see orthogonalError, same procedure except for complex numbers
addpath('~/software/gem/gem');
values_n = [2 4 8 16 32 64]; % use powers of two, so that nSamples is integer
errors = [];
nCoeffs = values_n(end)^2;
for i = 1:length(values_n)
    n = values_n(i);
    nSamples = nCoeffs / (n*n); % we want the same number of errors for each size
    err = [];
    fprintf('Computing for n=%d\n', n);
    for j = 1:nSamples
        % sample
        X = randn(n, n) + 1i*randn(n, n);
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
save -v7 unitaryError.mat errors values_n