classdef UnitaryGroup < replab.CompactGroup
% Describes the group of n x n unitary (complex) matrices
    
    properties
        n % integer: Dimension of the unitary group
    end
    
    properties (Access = protected)
        parent % replab.domain.Matrices: Domain of square complex matrices
    end
    
    methods
        
        function self = UnitaryGroup(n)
            self.n = n;
            self.parent = replab.domain.Matrices('C', n, n);
            self.identity = replab.eye_(n);
        end
        
        %% Str methods
        
        function s = headerStr(self)
            s = sprintf('%d x %d unitary matrices', self.n, self.n);
        end
        
        %% Domain methods
        
        function b = eqv(self, X, Y)
            b = self.parent.eqv(X, Y);
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
            X = self.parent.sample;
            [Q, R] = qr(X);
            R = diag(diag(R)./abs(diag(R)));
            X = Q*R;
        end

        %% Representation methods
        
        function rep = naturalRep(self)
        % Returns the natural representation of this unitary group
            rep = replab.Rep.lambda(self, 'C', self.n, true, @(u) u, @(u) u');
        end

    end

end
