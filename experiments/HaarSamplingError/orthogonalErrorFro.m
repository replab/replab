% Estimates the error on the coefficients of an orthogonal matrix sampled from the Haar measure
% using the technique in http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
%
% We sample a real matrix from the Haar measure, orthogonal up to double floating point precision
% We then compute the bound due to Felipe using GEM (to be sure)

addpath('~/software/gem/gem');
values_n = [2 4 8 16 32 64 128 256]; % use powers of two, so that nSamples is integer
nSamples = 64;
errors = zeros(nSamples, length(values_n));
for i = 1:length(values_n)
    n = values_n(i);
    fprintf('Computing for n=%d\n', n);
    for j = 1:nSamples
        % sample
        X = randn(n, n);
        [Q, R] = qr(X);
        F = diag(diag(R)./abs(diag(R)));
        Q = gem(Q*F);
        errors(j,i) = double(norm(Q*Q' - eye(n), 'fro'));
    end
end
save -v7 orthogonalErrorFro.mat errors values_n
