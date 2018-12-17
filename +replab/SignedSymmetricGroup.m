classdef SignedSymmetricGroup < replab.cat.Group
   
    properties
        n;
    end
    
    methods
        
        function self = SignedSymmetricGroup(n)
            self.n = n;
            self.parentOption = [];
            self.identity = 1:n;
        end
        
        function s = str(self)
            s = sprintf('Signed symmetric group acting on [-%d..1, 1..%d]', self.n, self.n);
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function h = hash(self, x)
            h = replab.cat.Domain.hashIntegers(x);
        end
        
        function s = sample(self)
            s = randperm(self.n) .* (randi([0 1], 1, self.n)*2-1);
        end
        
        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end
        
        function y = inverse(self, x)
            y = zeros(1, self.n);
            xAbs = abs(x);
            y(xAbs) = 1:self.n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end
        
        function A = naturalAction(self)
            n = self.n;
            desc = sprintf('Signed permutations acting on integers [-%d..-1, 1..%d]', n, n);
            P = replab.cat.Domain.integerRange(1, n);
            A = replab.cat.BSGSActionFun(desc, self, P, ...
                                         @(g, p) g(abs(p))*sign(p), ...
                                         @(g) replab.SignedPerm.findMovedPoint(g));

        end
        
        function phi = permIsomorphism(self)
            fun = @(x) replab.SignedPerm.toPerm(x);
            target = replab.SymmetricGroup(2*self.n);
            phi = replab.cat.GroupMorphismFun(fun, self, target);
        end
        
        function R = naturalRepresentation(self)
            fun = @(x) replab.SignedPerm.toMatrix(x);
            target = replab.SignedPermutationMatrices(self.n);
            R = replab.RepFun(fun, self, target);
        end
        
    end

end
