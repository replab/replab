classdef Isomorphism < replab.FiniteIsomorphism
% Describes the isomorphism from a group to an isomorphic abstract group

    methods

        function self = Isomorphism(source)
            self.source = source;
            self.target = replab.AbstractGroup(source.generatorNames, source.relatorsFlat, 'order', source.order, 'permutationGenerators', source.permutationGroup.generators);
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            t = self.source.factorizeWord(s);
        end

        % Isomorphism

        function s = preimageElement(self, t)
            s = self.source.imageWord(t);
        end

    end

end
