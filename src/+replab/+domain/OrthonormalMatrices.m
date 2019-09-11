classdef OrthonormalMatrices < replab.CompactGroup
% Describes the group of n x n orthonormal (real) matrices
    
    properties
        n; % size
    end

    properties (Access = protected)
        parent_; % real matrices
    end

    methods
        
        function self = OrthonormalMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.RealMatrices(n, n);
            if replab.Settings.useSparse
                self.identity = speye(n);
            else
                self.identity = eye(n);
            end
        end
        
        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d orthonormal matrices', self.n, self.n);
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.sampleUniformly;
        end

        %% Monoid methods
        
        function Z = compose(self, X, Y)
            Z = X * Y;
        end
        
        %% Group methods
        
        function XInv = inverse(self, X)
            XInv = X';
        end
        
        %% CompactGroup methods
        
        function X = sampleUniformly(self)
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            X = self.parent_.sample;
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

        %% Representations
        
        function rep = naturalRep(self)
        % Returns the natural representation of this orthonormal group
            rep = replab.Rep(self, 'R', self.n, @(o) o);
        end
        
    end

end
