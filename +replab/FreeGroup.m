classdef FreeGroup < replab.FinitelyGeneratedGroup & replab.cat.Group
    
    properties (SetAccess = protected)
        description;
        n; % number of generators
        canEqv = true;
        canHash = true;
        canSample = true;
    end
    
    methods
        
        function self = FreeGroup(n)
            self.description = sprintf('Free group on %d generators', self.n);
            self.n = n;
            generators = cell(1, n);
            for i = 1:n
                generators{i} = replab.Word.generator(i);
            end
            self.generators = generators;
            self.identity = replab.Word.identity;
        end
        
        % Domain methods
        
        function b = eqv(self, w1, w2)
            if length(w1.indices) ~= length(w2.indices)
                b = false;
                return
            end
            if length(w1.indices) == 0
                b = true;
                return
            end
            b = isequal(w1.indices, w2.indices) && isequal(w1.exponents, w2.exponents);
        end
        
        function h = hash(self, w)
            h = 31*replab.cat.Domain.hashIntegers(w.indices) + replab.cat.Domain.hashIntegers(w.exponents);
        end
        
        function W = sample(self, maxLength)
            if nargin < 2
                maxLength = 10;
            end
            l = randi(maxLength + 1) - 1;
            % exponents are +1 or -1
            e = randi(2, 1, l);
            e(e == 2) = -1;
            g = randi(self.nGenerators, 1, l); % generator indices
            W = replab.Word.fromIndicesAndExponents(g, e); % reduce
        end
        
        function W = compose(self, w1, w2)
            W = w1 * w2;
        end
        
        function W = inverse(self, w)
            W = inv(w);
        end
        
    end
    
end
