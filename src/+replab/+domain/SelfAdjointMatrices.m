classdef SelfAdjointMatrices < replab.domain.VectorSpace
% Describes the vector space of n x n symmetric/Hermitian matrices
    
    properties
        n % integer: Matrix size
    end
    
    properties (Access = protected)
        parent;
    end
    
    methods
        
        %% Own methods
        
        function self = SelfAdjointMatrices(field, n)
            self.n = n;
            self.parent = replab.domain.Matrices(field, n, n);
        end

        %% Str methods
        
        function s = headerStr(self)
            if self.overR
                s = sprintf('%d x %d symmetric real matrices', self.n, self.n);
            elseif self.overC
                s = sprintf('%d x %d Hermitian complex matrices', self.n, self.n);
            end
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end
        
        function X = sample(self)
            if self.overR
                % Generates a symmetric matrix with measure invariant under orthogonal transformations,
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
            elseif self.overC
                % Generates a Hermitian matrix with measure invariant under unitary transformations,
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

end
