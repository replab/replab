classdef PermMorphism < replab.Morphism

    properties (SetAccess = protected)
        chain % (`+replab.bsgs.ChainWithImages`): BSGS chain corresponding to the morphism
    end

    methods

        function self = PermMorphism(source, target, chain)
            assert(chain.n == source.domainSize);
            self.source = source;
            self.target = target;
            self.chain = chain;
        end

        function t = image(self, s)
            t = self.chain.image(s);
        end

    end

    methods (Static)

        function m = byImages(source, target, images)
            chain = replab.bsgs.ChainWithImages.make(source.domainSize, target, source.generators, images, [], source.chain.base, source.order);
            m = replab.mrp.PermMorphism(source, target, chain);
        end

    end

end
