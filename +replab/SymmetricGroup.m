classdef SymmetricGroup < replab.cat.Group
   
    properties
        domainSize;
    end
    
    methods
        
        function self = SymmetricGroup(domainSize)
            self.domainSize = domainSize;
            self.parentOption = [];
            self.identity = 1:domainSize;
        end
        
        function s = str(self)
            s = sprintf('Symmetric group acting on %d elements', self.domainSize);
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function h = hash(self, x)
            h = replab.cat.Domain.hashIntegers(x);
        end
        
        function s = sample(self)
            s = randperm(self.domainSize);
        end
        
        function z = compose(self, x, y)
            z = x(y);
        end
        
        function y = inverse(self, x)
            y = zeros(1, self.domainSize);
            y(x) = 1:self.domainSize;
        end
        
        function A = naturalAction(self)
            n = self.domainSize;
            desc = sprintf('Permutations of size %d acting on integers [1...%d]', n, n);
            P = replab.cat.Domain.integerRange(1, n);
            A = replab.cat.BSGSActionFun(desc, self, P, ...
                                         @(g, p) g(p), ...
                                         @(g) replab.Perm.findMovedPoint(g));
        end
        
        function R = naturalRepresentation(self)
            fun = @(x) replab.Perm.toMatrix(x);
            target = replab.PermutationMatrices(self.domainSize);
            R = replab.RepFun(fun, self, target);
        end
        
    end

end
