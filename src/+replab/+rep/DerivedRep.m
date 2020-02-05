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
            self.parent = parent;
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = parent.dimension;
            self.isUnitary = parent.isUnitary;
            self.conjugate = conjugate;
            self.inverse = inverse;
            self.transpose = transpose;
        end

        function s = headerStr(self)
            els = {};
            if self.conjugate
                els{1,end+1} = 'conjugate';
            end
            if self.inverse
                els{1,end+1} = 'inverse';
            end
            if self.transpose
                els{1,end+1} = 'transpose';
            end
            if isempty(els)
                els = 'Untransformed derived representation';
            else
                s = strjoin(els, ' ');
                s(1) = upper(s(1));
                s = [s ' derived representation'];
            end
        end

        % Rep

        function rho = image(self, g)
            if self.inverse
                img = self.parent.inverseImage(g);
            else
                img = self.parent.image(g);
            end
            if self.conjugate
                img = conj(img);
            end
            if self.transpose
                img = img.';
            end
        end

        function rho = inverseImage(self, g)
            if self.inverse
                img = self.parent.image(g);
            else
                img = self.parent.inverseImage(g);
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
