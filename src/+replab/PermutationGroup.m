classdef PermutationGroup < replab.NiceFiniteGroup
% A base class for all permutation groups

    properties (SetAccess = protected)
        domainSize % integer: The integer ``d``, as this group acts on ``{1, ..., d}``
    end

    properties (Access = protected)
        chain_ % (+replab.+bsgs.Chain): Stabilizer chain describing this permutation group
    end

    methods (Access = protected)

        function chain = computeChain(self)
            for i = 1:self.nGenerators
                assert(~self.isIdentity(self.generators{i}), 'Generators cannot contain the identity');
            end
            chain = replab.bsgs.Chain.make(self.domainSize, self.generators, [], self.order_);
        end

        function dec = computeDecomposition(self)
            c = self.chain;
            k = c.length;
            T = cell(1, k);
            for i = 1:k
                Ui = c.U{i};
                m = size(Ui, 2);
                Ti = cell(1, m);
                for j = 1:m
                    Ti{j} = Ui(:,j)';
                end
                T{i} = Ti;
            end
            dec = replab.FiniteGroupDecomposition(self, T);
        end

    end

    methods (Static)

        function pg = fromChain(chain, parent)
            if nargin < 2
                parent = [];
            end
            pg = replab.PermutationGroup(chain.n, chain.strongGenerators, chain.order, parent, chain);
        end

    end

    methods

        function self = PermutationGroup(domainSize, generators, order, parent, chain)
        % Constructs a permutation group
        %
        % Args:
        %   domainSize (integer): Size of the domain
        %   generators (cell(1,\*) of permutation): Group generators
        %   order (vpi, optional): Order of the group
        %   parent (`+replab.PermutationGroup`, optional): Parent of this group if known,
        %                                                  or ``'self'`` if this group is its own parent
        %   chain (`+replab.+bsgs.Chain`): BSGS chain describing the group
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
            self.generators = generators;
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            if nargin < 4
                parent = [];
            end
            if isempty(parent)
                self.parent = replab.S(domainSize);
            elseif isequal(parent, 'self')
                self.parent = self;
            else
                self.parent = parent;
            end
            if nargin > 4 && ~isempty(chain)
                self.chain_ = chain;
            end
            self.niceGroup_ = self;
            self.niceInverseMonomorphism_ = replab.Morphism.lambda(self, self, @(x) x);
        end

        function c = chain(self)
        % Returns the stabilizer chain corresponding to this permutation group
        %
        % Returns:
        %   `+replab.+bsgs.Chain`: Stabilizer chain
            if isempty(self.chain_)
                self.chain_ = self.computeChain;
            end
            c = self.chain_;
        end

        %% Domain methods

        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        function h = hash(self, x)
            h = replab.Domain.hashVector(x);
        end

        %% Monoid methods

        function z = compose(self, x, y)
            z = x(y);
        end

        %% Group methods

        function y = inverse(self, x)
            n = length(x);
            y = zeros(1, n);
            y(x) = 1:n;
        end

        function z = composeWithInverse(self, x, y)
            z = zeros(1, length(x));
            z(y) = x;
        end

        function z = leftConjugate(self, x, y)
            z = zeros(1, length(x));
            % x y xInv
            z(x) = x(y);
        end


        %% NiceFiniteGroup methods

        function res = hasSameParentAs(self, rhs)
            res = isa(rhs, 'replab.PermutationGroup') && (self.parent.domainSize == rhs.parent.domainSize);
        end

        function p = niceMonomorphismImage(self, p)
            p = p;
        end

        function grp = subgroup(self, generators, order)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators (row cell array): List of generators given as a permutations in a row cell array
        %   order (vpi, optional): Argument specifying the group order, if given can speed up computations
        %
        % Returns:
        %   +replab.PermutationGroup: The constructed permutation subgroup
            if nargin < 3
                order = [];
            end
            grp = replab.PermutationGroup(self.domainSize, generators, order, self.parent);
        end

        function o = elementOrder(self, p)
            orbits = replab.Partition.permutationsOrbits(p);
            orders = unique(orbits.blockSizes);
            o = 1;
            for i = 1:length(orders)
                o = lcm(o, orders(i));
            end
        end

        function res = isCyclic(self)
            if self.nGenerators <= 1
                res = true;
            elseif ~self.isCommutative
                res = false;
            else
                pds = factor(self.order);
                if length(pds) == 1 % group of prime order is cyclic
                    res = true;
                    return
                end
                pds = unique(pds);
                assert(all(pds <= 2^53-1)); % to be sure (otherwise can a BSGS have been computed?)
                pds = double(pds);
                for p = pds
                    newGens = cellfun(@(g) self.composeN(g, p), self.generators, 'uniform', 0);
                    newGens = newGens(~cellfun(@(g) self.isIdentity(g), newGens));
                    if self.subgroup(newGens).order*p ~= self.order
                        res = false;
                        return
                    end
                end
                res = true;
            end
        end


        function T = leftTransversal(self, subgroup)
            T = cellfun(@(g) self.inverse(g), self.rightTransversal(subgroup), 'uniform', 0);
        end

        function T = rightTransversal(self, subgroup)
            T = {};
            % equivalent to the non recursive code
            %          function rt = printRightTransversals(self, subgroup, g, subchain, l)
            %              if nargin == 2
            %                  g = self.identity;
            %                  subchain = subgroup.chain;
            %                  l = 1;
            %              end
            %              if l == self.chain.length + 1
            %                  g
            %              else
            %                  for b = self.chain.Delta{l}
            %                      bg = g(b);
            %                      orbit = subchain.orbitUnderG(1, bg);
            %                      if bg == min(orbit)
            %                          self.printRightTransversals(subgroup, g(self.chain.u(l, b)), subchain.stabilizer(bg), l + 1);
            %                      end
            %                  end
            %              end
            %          end
            n = self.domainSize;
            L = self.chain.length;
            stackG = zeros(n, L+1); % current group element
            stackG(:,1) = 1:n;
            stackS = cell(1, L+1);
            stackS{1} = subgroup.chain; % subgroup chain
            stackB = cell(1, L+1);
            stackB{1} = self.chain.Delta{1}; % remaining points in orbit to check
            l = 1;
            while l >= 1
                g = stackG(:,l)';
                B = stackB{l};
                while 1
                    B = stackB{l};
                    if isempty(B)
                        l = l - 1;
                        break
                    else
                        b = B(end);
                        stackB{l} = B(1:end-1);
                        bg = g(b);
                        orbit = stackS{l}.orbitUnderG(1, bg);
                        if bg == min(orbit)
                            if l == L
                                T{1,end+1} = g(self.chain.u(l, b));
                            else
                                stackS{l+1} = stackS{l}.stabilizer(bg);
                                stackG(:,l+1) = g(self.chain.u(l, b)); % compose transversal element
                                stackB{l+1} = self.chain.Delta{l+1};
                                l = l + 1;
                                break
                            end
                        end
                    end
                end
            end
        end

        function c = centre(self)
            c = self.centralizerGroup(self);
        end

        function res = union(self, other)
            assert(self.hasSameParentAs(other));
            res = self;
            for i = 1:other.nGenerators
                res = res.closure(other.generator(i));
            end
        end

        function res = intersection(self, other)
            assert(self.hasSameParentAs(other));
            if self.order > other.order
                res = other.intersection(self);
            else
                prop = @(g) other.contains(g);
                c = replab.bsgs.subgroupSearch(self.chain, prop);
                res = replab.PermutationGroup.fromChain(c);
            end
        end

        function c = centralizerElement(self, other)
            c = self.centralizerGroup(self.subgroup({other}));
        end

        function c = centralizerGroup(self, other)
            c = replab.bsgs.Centralizer(self, other).subgroup;
        end

        function sub = unorderedPartitionStabilizer(self, partition)
        % Computes the subgroup that leaves the given unordered partition invariant
        %
        % Example:
        %    >>> G = replab.S(4).unorderedPartitionStabilizer(replab.Partition.fromBlockIndices([1 1 2 2]));
        %    >>> G == replab.S(4).subgroup({[2 1 3 4] [3 4 1 2]})
        %        1
        %
        % Args:
        %   partition (`.Partition`): Unordered partition to leave invariant
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup that stabilizes the unordered partition
            n = partition.n;
            assert(n == self.domainSize);
            % sort the blocks from the smallest to the biggest
            blocks = partition.blocks;
            [~, I] = sort(cellfun(@length, blocks));
            blocks = blocks(I);
            % initialize the variables that will be captured by the function handles
            blockSizes = zeros(1, partition.n);
            blockIndex = zeros(1, partition.n);
            base = [];
            tests = cell(1, n);
            base = zeros(1, n);
            k = 1;
            for i = 1:length(blocks)
                blockSizes(blocks{i}) = length(blocks{i});
                blockIndex(blocks{i}) = i;
            end
            % create the tests
            for i = 1:length(blocks)
                block = blocks{i};
                s = length(block);
                b = block(1);
                base(k) = b;
                % the first element of a new block needs to be mapped to a block of the same size;
                % we pass on the index of that block we are mapped to as data
                tests{1,k} = @(g, data) deal(blockSizes(g(b)) == s, blockIndex(g(b)));
                k = k + 1;
                for j = 2:s
                    b = block(j);
                    base(k) = b;
                    % for the subsequent points in the same block, we verify that we still
                    % map to the same block
                    tests{1,k} = @(g, data) deal(blockIndex(g(b)) == data, data);
                    k = k + 1;
                end
            end
            c = self.chain.mutableCopy;
            c.baseChange(base);
            isConstant = @(x) all(x == x(1));
            prop = @(g) all(cellfun(@(b) isConstant(blockIndex(g(b))), blocks));
            subchain = replab.bsgs.subgroupSearch(c, prop, tests, []);
            sub = replab.PermutationGroup.fromChain(subchain, self.parent);
        end

        function sub = orderedPartitionStabilizer(self, partition)
        % Computes the subgroup that leaves the given ordered partition invariant
        %
        % The subgroup maps every block of the partition to itself under permutation of the
        % domain elements.
        %
        %
        % Example:
        %    >>> G = replab.S(4).orderedPartitionStabilizer(replab.Partition.fromBlockIndices([1 1 2 2]));
        %    >>> G == replab.S(4).subgroup({[2 1 3 4] [1 2 4 3]})
        %        1
        %
        % Args:
        %   partition (`.Partition`): Ordered partition to leave invariant
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup that stabilizes the ordered partition
            n = partition.n;
            assert(n == self.domainSize);
            blocks = partition.blocks;
            % sort the blocks from the smallest to the biggest
            [~, I] = sort(cellfun(@length, blocks));
            blocks = blocks(I);
            blockIndex = zeros(1, partition.n);
            base = [];
            tests = cell(1, n);
            base = zeros(1, n);
            k = 1;
            for i = 1:length(blocks)
                block = blocks{i};
                blockIndex(block) = i;
                for j = 1:length(block)
                    b = block(j);
                    base(k) = b;
                    % the tests express that every block is mapped to itself
                    tests{1,k} = @(g, data) deal(blockIndex(g(b)) == i, []);
                    k = k + 1;
                end
            end
            c = self.chain.mutableCopy;
            c.baseChange(base);
            prop = @(g) isequal(blockIndex(g), blockIndex);
            subchain = replab.bsgs.subgroupSearch(c, prop, tests, []);
            sub = replab.PermutationGroup.fromChain(subchain, self.parent);
        end

        function g = findPermutationTo(self, v, w, vStabilizer, wStabilizer)
        % Finds the permutation that relates two vectors, if it exists
        %
        % It returns the ``g`` such that ``w == v(g)``.
        %
        % Args:
        %   v (double(1,domainSize)): Row vector
        %   w (double(1,domainSize)): Another row vector
        %   vStabilizer (`.PermutationGroup`, optional): Vector stabilizer of ``v`` (or subgroup thereof)
        %   wStabilizer (`.PermutationGroup`, optional): Vector stabilizer of ``w`` (or subgroup thereof)
        %
        % Returns:
        %   permutation: The permutation ``g`` such that ``v == w(g)``
            if nargin < 4 || isequal(vStabilizer, [])
                vStabilizer = self.vectorStabilizer(v);
            end
            if nargin < 5 || isequal(wStabilizer, [])
                wStabilizer = self.vectorStabilizer(w);
            end
            chain = self.chain.mutableCopy;
            base = chain.base;
            tests = cell(1, length(base));
            for l = 1:length(base)
                tests{l} = @(g, data) deal(v(base(l)) == w(g(base(l))), []);
            end
            gInv = replab.bsgs.backtrackSearch(chain, @(x) isequal(v, w(x)), tests, [], vStabilizer.chain, wStabilizer.chain);
            g = self.inverse(gInv);
        end

        function sub = vectorStabilizer(self, vector)
        % Returns the permutation subgroup that leaves a given vector invariant
        %
        % Args:
        %   vector (double(1,domainSize)): Vector to stabilize under permutation
        %
        % Returns:
        %   `.PermutationGroup`: The subgroup of this group leaving ``vector`` invariant
            vector = vector(:).';
            v = unique(vector);
            c = arrayfun(@(x) sum(vector == x), v);
            [~, I] = sort(c);
            v = v(I);
            sub = self;
            for i = v
                sub = sub.setwiseStabilizer(find(vector == i));
            end
        end

        function s = setwiseStabilizer(self, set)
        % Returns the subgroup that stabilizes the given set as a set
        %
        % i.e. for this group ``G``, it returns ``H = {g \in G : g(set) = set}``
        %
        % Example:
        %   >>> G = replab.S(4).subgroup({[3 1 2 4] [1 4 2 3]});
        %   >>> H = G.setwiseStabilizer([1 2]);
        %   >>> H == replab.S(4).subgroup({[2 1 4 3]})
        %       1
        %
        % Example:
        %   >>> G = replab.S(10).subgroup({[1,3,2,10,9,8,6,5,7,4], [1,4,3,2,5,6,7,8,9,10]});
        %   >>> H = G.setwiseStabilizer([2 3]);
        %   >>> H.order
        %       10
        %
        % Args:
        %   set (integer(1,\*)): The set to stabilize
        %
        % Returns:
        %   `+replab.PermutationGroup`: The subgroup that stabilizes the set
            mask = false(1, self.domainSize);
            mask(set) = true;
            c = self.chain.mutableCopy;
            c.baseChange(set);
            tests = cell(1, length(set));
            for i = 1:length(tests)
                tests{i} = @(g, data) deal(mask(g(set(i))), []);
            end
            prop = @(g) all(mask(g(set)));
            subchain = replab.bsgs.subgroupSearch(c, prop, tests, []);
            s = replab.PermutationGroup.fromChain(subchain, self.parent);
        end

        function res = closureGroup(self, G)
            c = self.chain.mutableCopy;
            for i = 1:G.nGenerators
                if c.stripAndAddStrongGenerator(G.generator(i))
                    c.randomizedSchreierSims([]);
                end
            end
            c.makeImmutable;
            res = replab.PermutationGroup.fromChain(c, self.parent);
        end

        function res = closureElement(self, g)
            c = self.chain.mutableCopy;
            if c.stripAndAddStrongGenerator(g)
                c.randomizedSchreierSims([]);
            end
            c.makeImmutable;
            res = replab.PermutationGroup.fromChain(c, self.parent);
        end

        function nc = normalClosure(self, rhs)
            chain = replab.bsgs.Chain(self.domainSize);
            generators = {};
            toCheck = rhs.generators;
            while ~isempty(toCheck)
                rhsg = toCheck{end};
                toCheck = toCheck(1:end-1);
                for i = 1:self.nGenerators
                    gi = self.generator(i);
                    cm = self.leftConjugate(gi, rhsg);
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        toCheck{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            nc = replab.PermutationGroup(self.domainSize, generators, chain.order, self.parent, chain);
        end

        function sub = derivedSubgroup(self)
            nG = self.nGenerators;
            n = self.domainSize;
            chain = replab.bsgs.Chain(n);
            generators = {};
            for i = 1:nG
                gi = self.generator(i);
                for j = 1:nG
                    gj = self.generator(j);
                    cm = self.composeWithInverse(self.compose(gi, gj), self.compose(gj, gi));
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            % compute the normal closure
            toCheck = generators;
            while ~isempty(toCheck)
                h = toCheck{end};
                toCheck = toCheck(1:end-1);
                for i = 1:nG
                    gi = self.generator(i);
                    cm = self.leftConjugate(gi, h);
                    if chain.stripAndAddStrongGenerator(cm)
                        generators{1, end+1} = cm;
                        toCheck{1, end+1} = cm;
                        chain.randomizedSchreierSims([]);
                    end
                end
            end
            sub = replab.PermutationGroup(self.domainSize, generators, chain.order, self.parent, chain);
        end


        %% Methods specific to permutation groups

        function o = orbits(self)
        % Returns the partition of the domain into orbits under this group
        %
        % Permutation group orbits are also called domains of transitivity,
        % see https://www.encyclopediaofmath.org/index.php/Transitive_group
        %
        % Returns:
        %   `.Partition`: The orbit partition
            G = zeros(self.nGenerators, self.domainSize);
            for i = 1:self.nGenerators
                G(i, :) = self.generators{i};
            end
            o = replab.Partition.permutationsOrbits(G);
        end

        function sub = stabilizer(self, p)
        % Returns the subgroup that stabilizes a given point
            sub = replab.PermutationGroup.fromChain(self.chain.stabilizer(p), self.parent);
        end

        %% Group construction

        function w = wreathProduct(self, A)
        % Returns the wreath product of a compact group by this permutation group
        %
        % See https://en.wikipedia.org/wiki/Wreath_product
        %
        % Note that our notation is reversed compared to the Wikipedia page,
        % the permutation group is on the left hand side, as our convention
        % for semidirect product places the group acted upon on the right.
        %
        % Note that the return type depends on the argument type:
        % if ``A`` is a `.FiniteGroup`, the result will be a finite group too,
        % and if ``A`` is a `.NiceFiniteGroup`, the result will be of that type.
        %
        % Args:
        %   A (`.CompactGroup`): The group whose copies are acted upon
        %
        % Returns:
        %   `+replab.+wreathproduct.Common`: A wreath product group
            w = replab.wreathproduct.of(self, A);
        end

        %% Actions

        function A = naturalAction(self)
        % Returns the natural action of elements of this group on its domain
        %
        % This group natural domain is the set of integers ``{1..domainSize}``
        %
        % Returns:
        %   replab.Action: The natural action
            A = replab.perm.PermutationNaturalAction(self);
        end

        function A = vectorAction(self)
        % Returns the action of permutations on column vectors
        %
        % Acts on vectors of size `domainSize` by permuting their coefficients
        %
        % Returns:
        %   replab.Action: The vector action
            A = replab.perm.PermutationVectorAction(self);
        end

        function A = matrixAction(self)
        % Returns the simultaneous action of permutations on both rows and columns of square matrices
        %
        % Acts on matrices of size ``self.domainSize x self.domainSize``
        %
        % Returns:
        %   replab.Action: The matrix action
            A = replab.perm.PermutationMatrixAction(self);
        end

        function perm = indexRelabelingPermutation(self, g, indexRange)
        % Returns the permutation that acts by permuting tensor coefficients
        %
        % Let I = (i1, ..., id) be a sequence of indices, where d = self.domainSize
        % and 1 <= i1,...,id <= indexRange
        %
        % We enumerate elements of I by first incrementing id, then i_(d-1), etc...
        %
        % We compute the permutation of domain size ``indexRange^domainSize`` that acts on the
        % indices of I according to the argument ``g``.
        %
        % Args:
        %   g (permutation): Permutation of subindices
        %   indexRange (integer): Dimension of each subindex
        %
        % Returns:
        %   permutation: The permutation on the enumeration of indices
            n = self.domainSize;
            dims = indexRange * ones(1, n);
            perm = permute(reshape(1:prod(dims), dims), fliplr(n +  1 - g));
            perm = perm(:)';
        end

        function phi = indexRelabelingMorphism(self, indexRange)
        % Returns the morphism the permutation action of this group on tensor coefficients
        %
        % The tensor coefficients correspond to R^ir x R^ir ... (domainSize times)
        % where ir = indexRange
        %
        % See also:
        %   `+replab.PermutationGroup.indexRelabelingPermutation`
        %
        % Args:
        %   indexRange (integer): Dimension of each subindex
        %
        % Returns:
        %   function_handle: The permutation group homomorphism
            phi = @(g) self.indexRelabelingPermutation(g, indexRange);
        end

        %% Representation construction

        function rho = indexRelabelingRep(self, indexRange)
        % Representation that permutes the indices of a tensor
        %
        % It acts on the tensor space R^ir x R^ir ... (domainSize times)
        % where ir = indexRange, by permuting the indices.
        %
        % The representation returned is real.
        %
        % See also:
        %   `+replab.PermutationGroup.indexRelabelingPermutation`
        %
        % Args:
        %   indexRange (integer): Dimension of the tensor components/range of the subindices
        %
        % Returns:
        %   replab.Rep: The desired permutation representation
            rho = replab.rep.IndexRelabelingRep(self, indexRange);
        end

        function rho = naturalRep(self)
        % Returns the natural permutation representation of this permutation group
        %
        % Returns:
        %   replab.Rep: The (real) natural permutation representation
            rho = self.permutationRep(self.domainSize, self.generators);
        end

        function rho = standardRep(self)
        % Returns the standard representation of this permutation group
        %
        % It is an abuse of terminology as the "standard representation" is
        % the faithful $n-1$ dimensional representation of the symmetric group
        % acting on $n$ elements; but we can reuse that subrepresentation on
        % subgroups of the symmetric group.
        %
        % It corresponds to the representation orthogonal to the
        % trivial representation with basis ``[1, 1, ..., 1]'/sqrt(d)``
        %
        % Returns:
        %   `+replab.Rep`: The (real) standard representation
            [B_internal E_internal] = replab.sym.sageSpechtStandardBasis(self.domainSize);
            rho = self.naturalRep.subRep(B_internal, E_internal);
        end

        function rho = signRep(self)
        % Returns the sign representation of this permutation
            rho = replab.RepByImages.fromImageFunction(self, 'R', 1, @(g) replab.PermutationGroup.sign(g));
        end

    end

    methods(Static)

        function sign = sign(perm)
        % Returns the sign of a given permutation
        %
        % Args:
        %   perm (permutation): Vector representing a permutation (e.g. [3 2 1 4])
        %
        % Returns:
        %   integer: Sign of the permutation
            x = perm;
            n = length(x);
            oddOrEven = 0; %Records whether the total permutation is odd or even
            for i = 1:n
                if x(i) == 0 || x(i) == i %Skip over one cycles and numbers that have been cycled through
                    continue
                end
                cycleSize = -1; %The first element in a cycle isn't counted
                j = i;
                while x(j) ~= 0
                    pHold = x(j);
                    x(j) = 0;
                    j = pHold;
                    cycleSize = cycleSize + 1;
                end
                if cycleSize > 0
                    oddOrEven = oddOrEven + cycleSize; %At the end, this will match the parity (even/odd) of the permuation
                end
            end
            sign = (-1)^mod(round(oddOrEven),2); %Sign of permutation
            end
    end

end
