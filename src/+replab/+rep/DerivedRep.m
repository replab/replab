classdef DerivedRep < replab.Rep
% Representation derived by the means of complex conjugate, inverse and/or transpose

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation being transformed
        conjugate % (logical): Whether to take the complex conjugate
        inverse % (logical): Whether to take the inverse of the group element
        transpose % (logical): Whether to transpose the map
    end

    methods

        function self = DerivedRep(parent, conjugate, inverse, transpose)
            assert(inverse == transpose, ...
                   'Cannot use inverse & transpose independently, as our representations are left modules');
            % own properties
            self.parent = parent;
            self.conjugate = conjugate;
            self.inverse = inverse;
            self.transpose = transpose;
            % from replab.Rep, immutable
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = parent.dimension;
            % from replab.Rep, mutable
            self.isUnitary = parent.isUnitary;
            self.isTrivial = parent.isTrivial;
            self.isIrreducible = parent.isIrreducible;
        end

        function s = headerStr(self)
            s = self@replab.Rep.headerStr; % logic in parent class
        end

        % Rep

        function img = image_internal(self, g)
            if self.inverse
                img = self.parent.inverseImage_internal(g);
            else
                img = self.parent.image_internal(g);
            end
            if self.conjugate
                img = conj(img);
            end
            if self.transpose
                img = img.';
            end
        end

        function img = inverseImage_internal(self, g)
            if self.inverse
                img = self.parent.image_internal(g);
            else
                img = self.parent.inverseImage_internal(g);
            end
            if self.conjugate
                img = conj(img);
            end
            if self.transpose
                img = img.';
            end
        end

    end

end
