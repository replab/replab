classdef PermutationMatrices < replab.cat.Group
    
    properties
        n;
    end
    
    methods
        
        function self = PermutationMatrices(n)
            self.n = n;
            self.parentOption = [];
            self.identity = speye(n);
        end
        
        function s = str(self)
            s = sprintf('%d x %d permutation matrices', self.n, self.n);
        end
        
        function b = eqv(self, x, y)
            b = isequal(x, y);
        end
        
        function h = hash(self, x)
            h = replab.cat.Domain.hashIntegers(x);
        end
        
        function s = sample(self)
            P = randperm(self.n);
            s = sparse(P, 1:self.n, ones(1, self.n), self.n, self.n);
        end              
        
        function z = compose(self, x, y)
            z = x * y;
        end
        
        function y = inverse(self, x)
            y = inv(x);
        end
        
    end
    
end
