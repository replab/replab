classdef AutomorphismGroup < replab.NiceFiniteGroup

    properties (SetAccess = protected)
        object % (`.FiniteGroup`): Group

    end

    methods (Access = protected)
    end

    methods

        function self = AutomorphismGroup(object, generators)
            self.object = object;
            self.generators = generators;
            self.identity = replab.Automorphism.identity(object);
            self.type = self;
        end

        function s = objectSet(self)
            s = self.cached('objectSet', @() self.computeObjectSet);
        end

    end

    methods (Static)

        function A = innerAutomorphismGroup(object)
            generators = cellfun(@(g) replab.Automorphism.byConjugation(g, object), object.generators, 'uniform', 0);
            A = replab.AutomorphismGroup(object, generators);
        end

    end

    methods (Access = protected)

        function s = computeObjectSet(self)
            permGrp = self.object.niceMorphism.image;
            s = replab.perm.Set.fromPermutationGroup(permGrp);
            s.sort;
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(cellfun(@(g) self.object.eqv(x.imageElement(g), y.imageElement(g)), self.object.generators));
        end

        % Monoid

        function z = compose(self, x, y)
            zImg = cellfun(@(yi) x.imageElement(yi), y.generatorImages, 'uniform', 0);
            z = replab.Automorphism(self.object, zImg);
        end

        function z = inverse(self, x)
            z = x.inverse;
        end

        % NiceFiniteGroup

        function img = niceImage(self, g)
            A = g.niceAutomorphism;
            os = self.objectSet;
            n = os.nElements;
            img = zeros(1, n);
            for i = 1:n
                img(i) = os.find(A.imageElement(os.at(i)')');
            end
        end

    end

end
