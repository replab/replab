classdef UnitaryMatrices < replab.CompactGroup
% Describes the group of n x n unitary (complex) matrices
    
    properties
        n; % size
    end
    
    properties (Access = protected)
        parent_; % complex matrices
    end
    
    methods
        
        function self = UnitaryMatrices(n)
            self.n = n;
            self.parent_ = replab.domain.ComplexMatrices(n, n);
            if replab.Settings.useSparse
                self.identity = speye(n);
            else
                self.identity = eye(n);
            end
        end
        
        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d unitary matrices', self.n, self.n);
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.sampleUniformly;
        end
        
        % Semigroup/monoid/group
        
        function Z = compose(self, X, Y)
            Z = X * Y;
        end
        
        function XInv = inverse(self, X)
            XInv = X';
        end
        
        % Compact group
        
        function X = sampleUniformly(self)
        % see http://home.lu.lv/~sd20008/papers/essays/Random%20unitary%20[paper].pdf
            X = self.parent_.sample;
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

    end

end
