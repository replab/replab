classdef SignedPermutationGroup < replab.NiceFiniteGroup
% A base class for all signed permutation groups

    properties (SetAccess = protected)
        domainSize % d when this group acts on {-d..-1, 1..d}
    end

    methods

        function self = SignedPermutationGroup(domainSize, generators, order, type)
        % Constructs a signed permutation group
        %
        % Args:
        %   domainSize (integer): Size of the domain
        %   generators (cell(1,\*) of permutation): Group generators
        %   order (vpi, optional): Order of the group
        %   type (`+replab.SignedPermutationGroup`, optional): Type of this group if known,
        %                                                      or ``'self'`` if this group is its own type
            self.domainSize = domainSize;
            self.identity = 1:domainSize;
            self.generators = generators;
            if nargin > 2 && ~isempty(order)
                self.cache('order', order, '==');
            end
            if nargin < 4
                type = [];
            end
            if isempty(type)
                self.type = replab.SignedSymmetricGroup(domainSize);
            elseif isequal(type, 'self')
                self.type = self;
            else
                self.type = type;
            end
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.SignedPermutationGroupLaws(self);
        end

        % Domain

        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        % Monoid

        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end

        % Group

        function y = inverse(self, x)
            n = self.domainSize;
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end

        % NiceFiniteGroup methods

        function res = hasSameTypeAs(self, rhs)
            res = isa(rhs, 'replab.SignedPermutationGroup') && (self.type.domainSize == rhs.type.domainSize);
        end

        function p1 = nicePreimage(self, p)
            p1 = replab.SignedPermutation.fromPermutation(p);
        end

        function p1 = niceImage(self, p)
            p1 = replab.SignedPermutation.toPermutation(p);
        end

        function grp = niceSubgroup(self, generators, order, niceGroup)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators (row cell array): List of generators given as a permutations in a row cell array
        %   order (vpi, optional): Argument specifying the group order, if given can speed up computations
        %   niceGroup (`+replab.PermutationGroup`, optional): Image of this subgroup under the nice morphism
        %
        % Returns:
        %   +replab.SignedPermutationGroup: The constructed signed permutation subgroup
            if nargin < 4
                niceGroup = [];
            end
            if nargin < 3
                order = [];
            end
            grp = replab.SignedPermutationGroup(self.domainSize, generators, order, self.type);
            if ~isempty(niceGroup)
                grp.cache('niceGroup', niceGroup, '==');
            end
        end

        %% Methods specific to signed permutation groups

        function G = permutationPart(self)
        % Returns the permutation part of the current group
        %
        % Corresponds to the group image under the homomorphism `.elementPermutationPart`.
        %
        % Returns:
        %   replab.PermutationGroup: The corresponding permutation group
            newGenerators = cell(1, 0);
            for i = 1:self.nGenerators
                img = self.elementPermutationPart(self.generator(i));
                if ~self.isIdentity(img)
                    newGenerators{1, end+1} = img;
                end
            end
            G = replab.SymmetricGroup(self.domainSize).subgroup(newGenerators);
        end

        function p = elementPermutationPart(self, g)
        % Returns the permutation part of a signed permutation, by taking image absolute values
        %
        % Computes the permutation p that acts on 1...domainSize as p(i) = abs(g(i))
        %
        % Args:
        %   g (signed permutation): Signed permutation
        %
        % Returns:
        %   permutation: The permutation part of ``g``
            p = abs(g);
        end

        %% Actions

        function A = naturalAction(self)
        % Returns the action of elements of this group on {-d..-1 1..d}
        %
        % Here, d is self.domainSize
            A = replab.perm.SignedPermutationNaturalAction(self);
        end

        function A = vectorAction(self)
        % Returns the action of elements of this group on vectors
        %
        % Vectors are (self.domainSize)-dimensional vectors, and this permutes their coefficients and flips their signs
            A = replab.perm.SignedPermutationVectorAction(self);
        end

        function A = matrixAction(self)
        % Returns the action of elements of this group on d x d matrices
        %
        % Note that d = self.domainSize, and this acts by simultaneous permutations of their
        % rows and columns and flipping their signs
            A = replab.perm.SignedPermutationMatrixAction(self);
        end

        %% Representation construction

        function rho = naturalRep(self)
        % Natural representation on R^d of signed permutations on integers -d..-1, 1..d
            rho = self.signedPermutationRep(self.domainSize, self.generators);
        end

    end

    methods (Static)

        function G = trivial(n)
        % Returns
        end

        function G = of(varargin)
        % Constructs a nontrivial signed permutation group from the given generators
        %
        % If you do not know the number of generators in advance, and would like to handle the
        % case of a trivial group, use ``H = replab.SignedSymmetricGroup(n); H.subgroup(generators)`` instead.
        %
        % Example:
        %   >>> G = replab.SignedPermutationGroup.of([2 3 4 1], [4 3 2 1]);
        %   >>> G.order == 8
        %     1
        %
        % Args:
        %   varargin (cell(1,\*) of permutation): Group generators
        %
        % Returns:
        %   `+replab.SignedPermutationGroup`: The permutation group given as the closure of the generators
            assert(nargin > 0, 'Must be called with at least one generator');
            n = length(varargin{1});
            Sn = replab.SignedSymmetricGroup(n);
            G = Sn.subgroup(varargin);
        end

    end

end
