function M = sampleSymmetricMatrix(n)
% Generates a symmetric matrix with measure invariant 
% under orthogonal transformations,
% sampled from the Gaussian Orthogonal Ensemble, see
% http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
    M = zeros(n, n);
    for r = 1:n
        % diagonal elements are scaled up by sqrt(2)
        M(r, r) = randn * sqrt(2);
        % while other elements are standard normals
        M(r, r+1:end) = randn(1, n-r);
        M(r+1:end, r) = M(r, r+1:end);
    end
end
