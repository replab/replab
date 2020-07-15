classdef PermToPerm < replab.FiniteMorphism

    properties
        images % (cell(1,\*) of permutation): Images of the generators of the source group
    end

    methods (Access = protected)

        function z = join(self, x, y)
        % Computes the chain permutation corresponding to the permutations x, y
            n1 = length(x);
            n2 = length(y);
            z = zeros(1, n1+n2);
            z(1:n1) = x;
            z(n1+1:n1+n2) = y + n1;
        end

    end

    methods

        function self = PermToPerm(source, target, images)
            assert(isa(source, 'replab.PermutationGroup'));
            assert(isa(target, 'replab.PermutationGroup'));
            self.source = source;
            self.target = target;
            self.images = images;
        end

        function K = computeKernel(self)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            chain = self.inverseChain;
            l = find(chain.B <= n1, 1);
            if isempty(l)
                % the image group stabilizes
                K = self.source.trivialSubgroup;
            else
                base = chain.B(l:end);
                strongGens = chain.strongGeneratorsForLevel(l);
                orbitSizes = chain.orbitSizes;
                order = replab.util.multiplyIntegers(orbitSizes(l:end));
                generators = arrayfun(@(i) strongGens(1:n1,i)', 1:size(strongGens, 2), 'uniform', 0);
                chain = replab.bsgs.Chain.make(n1, generators, base, order);
                K = replab.PermutationGroup.fromChain(chain);
            end
        end

        function s = preimageRepresentative(self, t)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            chain = self.inverseChain;
            [h i] = self.chain.strip([1:n1 t+n1]);
            l = find(chain.B <= n1, 1);
            assert(i > l);
            sinv = h(1:n1);
            s(sinv) = 1:n1;
        end

        function S = preimageGroup(self, T)
            if T.isTrivial
                S = self.kernel;
            else
                preimages = cellfun(@(t) self.preimageRepresentative(t), T.generators, 'uniform', 0);
                S = self.kernel.closure(self.source.subgroup(preimages));
            end
        end

        function c = inverseChain(self)
            c = self.cached('inverseChain', @() self.computeInverseChain);
        end

        function c = computeInverseChain(self)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            c = self.chain.mutableCopy;
            c.baseChange(n1+1:n1+n2, true);
            c.makeImmutable;
        end

        function c = chain(self)
            c = self.cached('chain', @() self.computeChain);
        end

        function c = computeChain(self)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            S = arrayfun(@(i) self.join(self.source.generator(i), self.images{i}), 1:self.source.nGenerators, 'uniform', 0);
            chain = replab.bsgs.Chain.make(n1+n2, S, self.source.chain.base, self.source.order);
            chain = chain.mutableCopy;
            chain.baseChange(1:n1, true);
            chain.makeImmutable;
            c = chain;
        end

        function t = imageElement(self, s)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            chain = self.chain;
            [h i] = self.chain.strip([s n1+1:n1+n2]);
            assert(i == self.chain.length + 1);
            tinv = h(n1+1:n1+n2) - n1;
            t(tinv) = 1:n2;
        end

    end

end
