classdef DirectProductNiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        domainSizes % (integer(1,\*)): Size of the permutation images for each of the factors
        isomorphisms % (cell(1,\*) of `+replab.FiniteIsomorphism`): Isomorphisms into permutation groups for each factor group, that preserve the
    end

    methods

        function self = DirectProductNiceIsomorphism(sourceType, sourceGenerators, isomorphisms)
            self.isomorphisms = isomorphisms;
            domainSizes = cellfun(@(iso) iso.target.type.domainSize, isomorphisms);
            self.domainSizes = domainSizes;
            self.target = replab.S(sum(domainSizes));
        end

    end

    methods % Implementations

        function t = imageElement(self, s)
            shift = 0;
            p = [];
            n = length(self.isomorphisms);
            for i = 1:n
                pf = self.isomorphisms{i}.imageElement(s{i});
                p = [p pf+shift];
                shift = shift + length(pf);
            end
        end

        function s = preimageElement(self, t)
            shift = 0;
            n = length(self.isomorphisms);
            ds = self.domainSizes;
            s = cell(1, n);
            for i = 1:n
                pf = t((shift+1):shift+ds(i))-shift;
                shift = shift + ds(i);
                s{i} = self.isomorphisms{i}.preimageElement(pf);
            end
        end

    end

end
