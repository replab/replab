classdef DirectProductNiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        factors % (cell(1,\*) of `+replab.FiniteGroup`): Factors
        isomorphisms % (cell(1,\*) of `+replab.FiniteIsomorphism`): Order-preserving isomorphisms to permutation groups
        domainSizes % (integer(1,\*)): Size of the permutation images for each of the factors
    end

    methods

        function self = DirectProductNiceIsomorphism(factors)
            n = length(factors);
            assert(all(cellfun(@(x) isa(x, 'replab.FiniteGroup'), factors)));
            typeFactors = cellfun(@(g) g.type, factors, 'uniform', 0);
            type = replab.prods.DirectProductGroupType(typeFactors);
            generatorNames = cellfun(@(f) f.generatorNames, factors, 'uniform', 0);
            generatorNames = [generatorNames{:}];
            if length(unique(generatorNames)) ~= length(generatorNames) % names are not unique
                generatorNames = [];
            end
            isomorphisms = cellfun(@(f) f.orderPreservingPermutationIsomorphism, factors, 'uniform', 0);
            domainSizes = cellfun(@(iso) iso.target.type.domainSize, isomorphisms);
            sourceGenerators = cell(1, 0);
            targetGenerators = cell(1, 0);
            order = vpi(1);
            self.isomorphisms = isomorphisms;
            self.domainSizes = domainSizes;
            for i = 1:n
                iso = isomorphisms{i};
                G = factors{i};
                order = order * G.order;
                for j = 1:G.nGenerators
                    g = type.identity;
                    g{i} = G.generator(j);
                    sourceGenerators{1,end+1} = g;
                    targetGenerators{1,end+1} = self.imageElement(g);
                end
            end
            nice = replab.PermutationGroup(sum(domainSizes), targetGenerators, 'order', order, 'generatorNames', generatorNames);
            self.source = replab.prods.DirectProductGroup_finite(factors, type, sourceGenerators, nice, self);
            self.target = nice;
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            shift = 0;
            t = [];
            n = length(self.isomorphisms);
            for i = 1:n
                pf = self.isomorphisms{i}.imageElement(s{i});
                t = [t pf+shift];
                shift = shift + length(pf);
            end
        end

        % Isomorphism

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

        % NiceIsomorphism

        function l = sourceContains(self, s)
            l = false;
            for i = 1:length(self.factors)
                if ~self.factors{i}.contains(s{i})
                    return
                end
            end
            l = true;
        end

    end

end
