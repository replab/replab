function M = sampleHermitianMatrix(n)
% Generates a n x n Hermitian matrix with measure invariant
% under unitary transformations,
% sampled from the Gaussian Unitary Ensemble, see
% http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
    M = zeros(n, n);
    for r = 1:n
        M(r, r) = randn;
        M(r, r+1:end) = (randn(1, n-r) + randn(1, n-r)*1i)/sqrt(2);
        M(r+1:end, r) = conj(M(r, r+1:end));
    end    
end
