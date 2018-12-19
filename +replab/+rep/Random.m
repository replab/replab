classdef Random
    
    methods (Static)
        
        function M = realGaussian(d)
            M = randn(d);
        end
        
        function M = complexGaussian(d)
            M = (randn(d) + randn(d) * 1i)/sqrt(2);
        end
        
        function M = hermitianGaussian(d)
        % Generates a Hermitian matrix with measure invariant under unitary transformations,
        % sampled from the Gaussian Unitary Ensemble
        % see http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
        %
        % d is the dimension            
            M = zeros(d, d);
            for r = 1:d
                M(r, r) = randn;
                M(r, r+1:end) = (randn(1, d-r) + randn(1, d-r)*1i)/sqrt(2);
                M(r+1:end, r) = conj(M(r, r+1:end));
            end
        end
        
        function M = symmetricGaussian(d)
        % Generates a symmetric matrix with measure invariant under orthogonal transformations,
        % sampled from the Gaussian Orthogonal Ensemble
        % see http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
        %
        % d is the dimension
            M = zeros(d, d);
            for r = 1:d
                % diagonal elements are scaled up by sqrt(2)
                M(r, r) = randn * sqrt(2);
                % while other elements are standard normals
                M(r, r+1:end) = randn(1, d-r);
                M(r+1:end, r) = M(r, r+1:end);
            end
        end
                        
    end
    
end
