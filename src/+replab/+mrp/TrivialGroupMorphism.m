classdef TrivialGroupMorphism < replab.Morphism

    methods

        function self = CompositionMorphism(source, target)
            self.source = source;
            self.target = target;
        end

        function t = image(s)
            t = self.target.identity;
        end

    end

end
