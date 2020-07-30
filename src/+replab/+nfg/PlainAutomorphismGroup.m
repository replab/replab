classdef PlainAutomorphismGroup < replab.Group

    properties (SetAccess = protected)
        object % (`.FiniteGroup`): Group
    end

    methods

        function self = PlainAutomorphismGroup(object)
            self.object = object;
            self.identity = replab.nfg.PlainAutomorphism.identity(object);
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = all(cellfun(@(g) self.object.eqv(x.imageElement(g), y.imageElement(g)), self.object.generators));
        end

        function s = sample(self)
            s = replab.nfg.PlainAutomorphism.byConjugation(self.object.sample, self.object);
        end

        % Monoid

        function z = compose(self, x, y)
            zImg = cellfun(@(yi) x.imageElement(yi), y.imageSourceGenerators, 'uniform', 0);
            z = replab.nfg.PlainAutomorphism(self.object, zImg);
        end

        function z = inverse(self, x)
            zImg = cellfun(@(g) x.preimageRepresentative(g), x.source.generators, 'uniform', 0);
            z = replab.nfg.PlainAutomorphism(self.object, zImg);
        end

    end

end
