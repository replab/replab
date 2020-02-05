classdef SimilarRep < replab.Rep
% Describes a representation similar to a given representation under a change of basis
%
% We use left action convention, which means that
% ``image(g) = A * parent.image(g) * inv(A)``

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation this is derived from
        A % (double(*,*)): Change of basis matrix
        Ainv % (double(*,*)): Inverse of change of basis matrix
    end

    methods

        function self = SimilarRep(parent, A, Ainv)
            switch parent.field
              case 'R'
                assert(isreal(A) && isreal(Ainv), 'A real Rep can only be conjugated by a real orthonormal matrix');
              case 'C'
                assert(isa(A, 'double') && isa(Ainv, 'double'), 'A complex Rep can only be conjugated by a complex/real unitary matrix');
              otherwise
            end
            assert(size(A, 1) == parent.dimension);
            assert(size(A, 2) == parent.dimension);
            assert(isequal(size(A), size(Ainv)));
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = parent.dimension;
            self.A = A;
            self.Ainv = Ainv;
            self.parent = parent;
            self.irrepInfo = [];
        end

        function s = headerStr(self)
            s = 'Similar representation';
        end

        function rho = image(self, g)
            rho = self.A * self.parent.image(g) * self.Ainv;
        end

        function rho = inverseImage(self, g)
            rho = self.A * self.parent.inverseImage(g) * self.Ainv;
        end

    end

end
