classdef FiniteGroup < replab.CompactGroup & replab.GroupWithGenerators
% Describes a group with a finite number of elements
%
% More computational primitives are available when a nice monomorphism is provided, see `.NiceFiniteGroup`.
%
% This computational structure must provide:
%
% * An indexed family of the group elements that supports element seeking and retrieval
%
% * A decomposition of the finite group in a product of sets

    properties (Access = protected)
        randomBag_ % `.RandomBag`: Bag of elements used for random sampling
        order_ % vpi: Cached group order
        elements_ % `.IndexedFamily`: Cached indexed family of group elements
        decomposition_ % `.FiniteGroupDecomposition`: Cached decomposition of this group
    end

    methods (Access = protected) % Abstract

        function o = computeOrder(self)
        % Computes the result cached by self.order
            error('Abstract');
        end

        function E = computeElements(self)
        % Computes the result cached by self.elements
            error('Abstract');
        end

        function D = computeDecomposition(self)
        % Computes the result cached by self.decomposition
            error('Abstract');
        end

    end

    methods % Abstract (and cached) methods

        function o = order(self)
        % Returns the group order
        %
        % Returns:
        %   vpi: The group order
            if isempty(self.order_)
                self.order_ = self.computeOrder;
            end
            o = self.order_;
        end

        function E = elements(self)
        % Returns an indexed family of the group elements
        %
        % The order of elements in the family is not guaranteed due to the use of nondeterministic algorithms.
        %
        % Returns:
        %   `.IndexedFamily`: A space-efficient enumeration of the group elements
            if isempty(self.elements_)
                self.elements_ = self.computeElements;
            end
            E = self.elements_;
        end

        function D = decomposition(self)
        % Returns a decomposition of this group as a product of sets
        %
        % Returns:
        %   `.FiniteGroupDecomposition`: The group decomposition
            if isempty(self.decomposition_)
                self.decomposition_ = self.computeDecomposition;
            end
            D = self.decomposition_;
        end

    end

    methods

        %% Own methods

        function R = randomBag(self)
        % Returns an instance of the product-replacement algorithm data structure
        %
        % Returns:
        %   `.RandomBag`: The created random bag
            if isempty(self.randomBag_)
                self.randomBag_ = replab.fg.RandomBag(self, self.generators);
            end
            R = self.randomBag_;
        end

        function res = isCommutative(self)
        % Returns whether this group is commutative
            for i = 1:self.nGenerators
                gi = self.generator(i);
                for j = 1:i-1
                    gj = self.generator(j);
                    if ~self.eqv(self.compose(gi, gj), self.compose(gj, gi))
                        res = false;
                        return
                    end
                end
            end
            res = true;
            return
        end

        %% Str methods

        function names = hiddenFields(self)
            names = hiddenFields@replab.Group(self);
            names{1, end+1} = 'generators';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Group(self);
            for i = 1:self.nGenerators
                names{1, end+1} = sprintf('generator(%d)', i);
                values{1, end+1} = self.generator(i);
            end
        end

        %% Domain methods

        function g = sample(self)
            g = self.randomBag.sample;
        end

        %% CompactGroup methods

        function g = sampleUniformly(self)
            g = self.elements.sample;
        end

        %% Representations

        function rep = regularRep(self)
        % Returns the left regular representation of this group
        %
        % Warning:
        %   The choice of a basis for the regular representation is not deterministic,
        %   and results can vary between runs.
        %
        % Returns:
        %   `.Rep`: The left regular representation as a real permutation representation
            o = self.order;
            assert(o < 1e6);
            o = double(o);
            perms = cell(1, self.nGenerators);
            E = self.elements;
            for i = 1:self.nGenerators
                g = self.generator(i);
                img = zeros(1, o);
                for j = 1:o
                    img(j) = double(E.find(self.compose(g, E.at(j))));
                end
                perms{i} = img;
            end
            rep = self.permutationRep(o, perms);
        end

    end

end
