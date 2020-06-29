classdef PermPermMorphism < replab.Morphism

    properties (SetAccess = protected)
        chain % (`+replab.bsgs.Chain`): BSGS chain corresponding to the morphism
    end

    methods

        function self = PermPermMorphism(source, target, chain)
            assert(chain.n == source.domainSize + target.domainSize);
            self.source = source;
            self.target = target;
            self.chain = chain;
        end

        function t = image(self, s)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            [h i] = self.chain.strip([s n1+1:n1+n2]);
            assert(i == self.chain.length + 1);
            el = h(n1+1:n1+n2) - n1;
            t(1:n2) = el;
        end

    end

    methods (Static)

        function m = fromImages(source, images)
            nG = source.nGenerators;
            assert(nG > 1);
            n1 = source.domainSize;
            n2 = length(images{1});
            S = zeros(n1+n2, nG);
            for i = 1:nG
                S(1:n1, i) = source.generator(i);
                S(n1+1:n1+n2, i) = images{i} + n1;
            end
            chain = replab.bsgs.Chain.make(n1+n2, S, 1:n1);
            m = replab.mrp.PermPermMorphism(source, replab.S(n2), chain);
        end

    end

end
