classdef SimilarRep < replab.Rep
% Describes a representation similar to a given representation under a change of basis
%
% We use left action convention, which means that
% ``image(g) = A * parent.image(g) * inv(A)``

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation this is derived from
        A_internal % (double(*,*), may be sparse): Change of basis matrix
        Ainv_internal % (double(*,*), may be sparse): Inverse of change of basis matrix
    end

    methods

        function self = SimilarRep(parent, A_internal, Ainv_internal)
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
            self.A_internal = A_internal;
            self.Ainv_internal = Ainv_internal;
            self.parent = parent;
            self.irrepInfo = [];
        end

        function s = headerStr(self)
            s = 'Similar representation';
        end

        function rho = image_internal(self, g)
            rho = self.A_internal * self.parent.image_internal(g) * self.Ainv_internal;
        end

        function rho = inverseImage_internal(self, g)
            rho = self.A_inverse * self.parent.inverseImage_internal(g) * self.Ainv_internal;
        end

        function mat = A(self)
            mat = full(self.A_internal);
        end

        function mat = Ainv(self)
            mat = full(self.Ainv_internal);
        end

    end

end
