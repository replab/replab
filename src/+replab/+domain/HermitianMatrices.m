classdef HermitianMatrices < replab.Domain
% Describes the vector space of n x n (complex) Hermitian matrices
    
    properties
        n; % matrix size
    end
    
    properties (Access = protected)
        parent_;
    end
    
    methods
        
        function self = HermitianMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.ComplexMatrices(n, n);
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
        % Generates a n x n Hermitian matrix with measure invariant
        % under unitary transformations,
        % sampled from the Gaussian Unitary Ensemble, see
        % http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
            n = self.n;
            X = zeros(n, n);
            for r = 1:n
                X(r, r) = randn;
                X(r, r+1:end) = (randn(1, n-r) + randn(1, n-r)*1i)/sqrt(2);
                X(r+1:end, r) = conj(X(r, r+1:end));
            end
        end    
        
    end

end
