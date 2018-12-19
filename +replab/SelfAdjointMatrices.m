classdef SelfAdjointMatrices < replab.cat.Domain
    
    properties (SetAccess = protected)
        n;
        field;
        identity;
        canEqv = false;
        canHash = false;
        canSample = true;
    end
    
    methods
        
        function self = SelfAdjointMatrices(n, field)
            self.n = n;
            self.field = field;
            switch self.field
              case {'R15', 'C15'}
                self.identity = eye(n);
              case {'R7', 'C7'}
                self.origin = eye(n, 'single');
              otherwise
                error(sprintf('Unknown field %s', field));
            end
        end
        
        function s = str(self)
            s = sprintf('%d x %d dimensional self-adjoint in ', self.n, self.n, self.field);
        end
        
        function b = eqv(self, x, y)
            error('Cannot test floating point matrices for equality');
        end
        
        function h = hash(self, x)
            error('Cannot hash floating point matrices');
        end
        
        function s = sample(self)
        % TODO: make faster and corresponding to qdimsum GOE/GUE
            n = self.n;
            switch self.field
              case 'R7'
                s = single(randn(n, n));
              case 'C7'
                s = single(randn(n, n)) + 1i*single(randn(n, n));
              case 'R15'
                s = randn(n, n);
              case 'C15'
                s = randn(n, n) + 1i*randn(n, n);
            end
            s = (s + s')/2;
        end
        
    end
    
end
