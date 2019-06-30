classdef QuaternionHermitianMatrices < replab.Domain
% Describes the vector space of n x n quaternion Hermitian matrices
    
    properties
        n; % matrix size
    end
    
    properties (Access = protected)
        parent_;
    end
    
    methods
        
        function self = QuaternionHermitianMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.QuaternionMatrices(n, n);
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d quaternion Hermitian matrices', self.n, self.n);
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
            X = zerosq(n, n);
            for r = 1:n
                X(r, r) = randn;
                s1 = randn(1, n-r);
                si = randn(1, n-r);
                sj = randn(1, n-r);
                sk = randn(1, n-r);
                if r < n
                    X(r, r+1:end) = (s1*q1 + si*qi + sj*qj + sk*qk)/2;
                    X(r+1:end, r) = conj(X(r, r+1:end));
                end
            end
        end    
        
    end

end
