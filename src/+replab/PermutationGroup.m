classdef PermutationGroup < replab.NiceFiniteGroup
% A base class for all permutation groups

    properties (SetAccess = protected)
        domainSize % integer: The integer ``d``, as this group acts on ``{1, ..., d}``
    end

    methods

        function self = PermutationGroup(domainSize, generators, order, parent)
        % Constructs a permutation group
        %
        % Args:
        %   domainSize (integer): Size of the domain
        %   generators (cell(1,\*) of permutation): Group generators
        %   order (vpi, optional): Order of the group
        %   parent (`replab.PermutationGroup`, optional): Parent of this group if known,
        %                                                 or ``[]`` if this group is its own parent
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
            self.generators = generators;
            if nargin > 2 && ~isempty(order)
                self.order_ = order;
            end
            if nargin > 3
                if isempty(parent)
                    self.parent = self;
                else
                    self.parent = parent;
                end
            else
                self.parent = replab.S(domainSize);
            end
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
            n = self.domainSize;
            y = zeros(1, n);
            y(x) = 1:n;
        end

        %% NiceFiniteGroup methods

        function res = sameParentAs(self, rhs)
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
