classdef ConjugateRep < replab.Rep
% Conjugated representation U * rep * U'
%
% We use left action convention, which means that
% image(g) = U * parent.image(g) * U'

    properties (SetAccess = protected)
        U % unitary conjugation matrix
        parent % representation being conjugated
    end

    methods

        function self = ConjugateRep(U, parent)
            switch parent.field
              case 'R'
                assert(isreal(U), 'A real Rep can only be conjugated by a real orthonormal matrix');
              case 'C'
                assert(isa(U, 'double'), 'A complex Rep can only be conjugated by a complex/real unitary matrix');
              otherwise
            end
            assert(size(U, 1) == parent.dimension);
            assert(size(U, 2) == parent.dimension);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = parent.dimension;
            self.U = U;
            self.parent = parent;
            self.irrepInfo = [];
        end

        function s = headerStr(self)
            s = 'Conjugated representation';
        end

        function rho = image(self, g)
            rho = self.U * self.parent.image(g) * self.U';
        end

    end

end
