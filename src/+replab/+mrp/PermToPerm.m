classdef PermToPerm < replab.FiniteMorphism
% Describes a morphism between two permutation groups
%
% We use the approach in Section 4.5.3 of Derek Holt, Handbook of CGT

    properties (SetAccess = protected)
        preimages % (cell(1,\*) of permutation): Preimages in the source group
        images % (cell(1,\*) of permutation): Images in the target group
        imageElementFun % (function_handle or ``[]``): Optional function_handle that implements `.imageElement`
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

    methods (Access = protected)

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

        function c = computeInverseChain(self)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            c = self.chain.mutableCopy;
            c.baseChange(n1+1:n1+n2, true);
            c.makeImmutable;
        end

        function c = computeChain(self)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            S = arrayfun(@(i) self.join(self.preimages{i}, self.images{i}), 1:length(self.preimages), 'uniform', 0);
            chain = replab.bsgs.Chain.make(n1+n2, S, self.source.chain.base, self.source.order);
            chain = chain.mutableCopy;
            chain.baseChange(1:n1, true);
            chain.makeImmutable;
            c = chain;
        end

    end

    methods

        function self = PermToPerm(source, target, preimages, images, imageElementFun)
            assert(isa(source, 'replab.PermutationGroup'));
            assert(isa(target, 'replab.PermutationGroup'));
            self.source = source;
            self.target = target;
            self.preimages = preimages;
            self.images = images;
            self.imageElementFun = imageElementFun;
        end

        function c = inverseChain(self)
            c = self.cached('inverseChain', @() self.computeInverseChain);
        end

        function c = chain(self)
            c = self.cached('chain', @() self.computeChain);
        end

    end

    methods % Implementations

        function s = preimageRepresentative(self, t)
            n1 = self.source.domainSize;
            n2 = self.target.domainSize;
            el = [1:n1 t+n1];
            [h i] = self.inverseChain.strip(el);
            l = find([self.inverseChain.B n1] <= n1, 1);
            assert(i >= l);
            sinv = h(1:n1);
            s(sinv) = 1:n1;
        end

        function t = imageElement(self, s)
            if ~isempty(self.imageElementFun)
                f = self.imageElementFun;
                t = f(s);
                return
            end
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
