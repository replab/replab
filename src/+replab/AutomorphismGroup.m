classdef AutomorphismGroup < replab.Group
% Describes the automorphism group of a group

    properties (SetAccess = protected)
        object % (`.Group`): Group
    end

    methods

        function self = AutomorphismGroup(object)
            self.object = object;
            if isa(object, 'replab.FiniteGroup')
                self.identity = replab.FiniteIsomorphism.identity(object);
            else
                self.identity = replab.Isomorphism.identity(object);
            end
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(cellfun(@(g) self.object.eqv(x.imageElement(g), y.imageElement(g)), self.object.generators));
        end

        function s = sample(self)
            g = self.object.sample;
            s = self.object.innerAutomorphism(g);
        end

        % Monoid

        function z = compose(self, x, y)
            if isa(self.object, 'replab.FiniteGroup')
                zImg = cellfun(@(g) x.imageElement(y.imageElement(g)), self.object.generators, 'uniform', 0);
                z = self.object.isomorphismByImages(self.object, 'preimages', self.object.generators, 'images', zImg);
            else
                z = x.compose(y);
            end
        end

        function z = inverse(self, x)
            if isa(self.object, 'replab.FiniteGroup')
                zImg = cellfun(@(g) x.preimageElement(g), self.object.generators, 'uniform', 0);
                z = self.object.isomorphismByImages(self.object, 'preimages', self.object.generators, 'images', zImg);
            else
                z = x.inverse;
            end
        end

    end

end
