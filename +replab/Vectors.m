classdef Vectors < replab.cat.Domain
    
    properties
        n;
        field;
        origin;
    end
    
    methods
        
        function self = Vectors(n, field)
            self.n = n;
            self.field = field;
            switch self.field
              case {'R15', 'C15'}
                self.origin = zeros(n, 1);
              case {'R7', 'C7'}
                self.origin = zeros(n, 1, 'single');
              otherwise
                error(sprintf('Unknown field %s', field));
            end
            self.parentOption = [];
        end
        
        function s = str(self)
            s = sprintf('%d dimensional vectors in ', self.n, self.field);
        end
        
        function b = eqv(self, x, y)
            b = norm(x - y) < replab.Settings.eigTol;
        end
        
        function h = hash(self, x)
            error('Cannot hash floating point vectors');
        end
        
        function s = sample(self)
            n = self.n;
            switch self.field
              case 'R7'
                s = single(randn(n, 1));
              case 'C7'
                s = single(randn(n, 1)) + 1i*single(randn(n, 1));
              case 'R15'
                s = randn(n, 1);
              case 'C15'
                s = randn(n, 1) + 1i*randn(n, 1);
            end
        end
        
    end
    
end
