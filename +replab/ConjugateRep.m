classdef ConjugateRep < replab.Rep
% Conjugated representation A * rep * A'
%
% We use left action convention, which means that
% image(g) = A * parent.image(g) * A'
    
    properties (SetAccess = protected)
        A; % unitary conjugation matrix
        parent; % representation being conjugated
    end
    
    methods
                
        function self = ConjugateRep(A, parent)
            switch parent.field
              case 'R'
                assert(isreal(A), 'A real Rep can only be conjugated by a real orthonormal matrix');
              case 'C'
                assert(isa(A, 'double'), 'A complex Rep can only be conjugated by a complex/real unitary matrix');
              otherwise
            end
            assert(size(A, 1) == parent.dimension);
            assert(size(A, 2) == parent.dimension);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = parent.dimension;
            self.A = A;
            self.parent = parent;
        end

        function s = headerStr(self)
            s = 'Conjugated representation';
        end

        function rho = image(self, g)
            rho = self.A * self.parent.image(g) * self.A';
        end

        % TODO: optimize other methods
    end
    
end
