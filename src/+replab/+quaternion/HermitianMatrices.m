classdef HermitianMatrices < replab.Domain
% Describes the vector space of n x n quaternion Hermitian matrices
    
    properties
        n % integer: Matrix size
    end
    
    properties (Access = protected)
        parent;
    end
    
    methods
        
        function self = HermitianMatrices(n)
            self.n = n;
            self.parent = replab.quaternion.Matrices(n, n);
        end
        
        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d quaternion Hermitian matrices', self.n, self.n);
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
        end
        
        function X = sample(self)
        % Generates a n x n quaternion Hermitian matrix
        %
        % We adapt the procedure used for complex Hermitian
        % matrices but haven't proved that it is the correct
        % approach (factors could be off)
            n = self.n;
            X = replab.quaternion.H.zeros(n, n);
            for r = 1:n
                X(r, r) = randn;
                s1 = randn(1, n-r);
                si = randn(1, n-r);
                sj = randn(1, n-r);
                sk = randn(1, n-r);
                if r < n
                    X(r, r+1:end) = replab.quaternion.H(s1, si, sj, sk)/2;
                    X(r+1:end, r) = conj(X(r, r+1:end));
                end
            end
        end    

    end

end
