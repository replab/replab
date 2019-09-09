classdef SymmetricMatrices < replab.Domain
% Describes the vector space of n x n (real) symmetric matrices
    
    properties
        n; % matrix size
    end
    
    properties (Access = protected)
        parent_;
    end
    
    methods
        
        function self = SymmetricMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.RealMatrices(n, n);
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d Hermitian matrices', self.n, self.n);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
        % Generates a symmetric matrix with measure invariant 
        % under orthogonal transformations,
        % sampled from the Gaussian Orthogonal Ensemble, see
        % http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
            X = zeros(n, n);
            for r = 1:n
                % diagonal elements are scaled up by sqrt(2)
                X(r, r) = randn * sqrt(2);
                % while other elements are standard normals
                X(r, r+1:end) = randn(1, n-r);
                X(r+1:end, r) = X(r, r+1:end);
            end
        end
        
    end

end
