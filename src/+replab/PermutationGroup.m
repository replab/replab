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
                pds = unique(factor(self.order));
                assert(all(pds <= 2^53-1)); % to be sure, but unlikely (otherwise can a BSGS be computed?)
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


        function T = leftTransversals(self, subgroup)
            T = cellfun(@(g) self.inverse(g), self.rightTransversals(subgroup), 'uniform', 0);
        end

        function T = rightTransversals(self, subgroup)
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

        function c = centralizerElement(self, other)
        % Returns the centralizer of a group in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerElement([2 3 1 4]);
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Args:
        %   other (`+replab.PermutationGroup`): Permutation group with the same
            c = self.centralizerGroup(self.subgroup({other}));
        end

        function c = centralizerGroup(self, other)
        % Returns the centralizer of a group in this group
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerGroup(G.subgroup({[2 3 1 4]}));
        %   >>> C == replab.S(4).subgroup({[2 3 1 4]})
        %     1
        %
        % Example:
        %   >>> G = replab.S(4);
        %   >>> C = G.centralizerGroup(G.subgroup({[2 3 1 4] [2 1 3 4]}));
        %   >>> C == replab.S(4).trivialGroup
        %     1
        %
        % Args:
        %   other (`+replab.PermutationGroup`): Permutation group with the same
            c = replab.bsgs.Centralizer(self, other).subgroup;
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
                tests{i} = @(g) mask(g(set(i)));
            end
            prop = @(g) all(mask(g(set)));
            subchain = replab.bsgs.subgroupSearch(c, prop, set, tests);
            s = replab.PermutationGroup.fromChain(subchain, self.parent);
        end

        function res = closure(self, g)
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
        % trivial representation with basis [1, 1, ..., 1]'/sqrt(d)
        %
        % Returns:
        %   `+replab.Rep`: The (real) standard representation
            [B_internal E_internal] = replab.sym.sageSpechtStandardBasis(self.domainSize);
            rho = self.naturalRep.subRep(B_internal, E_internal);
        end

    end

end
