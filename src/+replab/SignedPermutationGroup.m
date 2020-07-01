classdef SignedPermutationGroup < replab.NiceFiniteGroup
% A base class for all signed permutation groups

    properties (SetAccess = protected)
        domainSize % d when this group acts on {-d..-1, 1..d}
    end

    methods

        function self = SignedPermutationGroup(domainSize, generators, order, parent)
        % Constructs a signed permutation group
        %
        % Args:
        %   domainSize (integer): Size of the domain
        %   generators (cell(1,\*) of permutation): Group generators
        %   order (vpi, optional): Order of the group
        %   parent (`replab.SignedPermutationGroup`, optional): Parent of this group if known,
        %                                                       or ``[]`` if this group is its own parent
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
                self.parent = replab.SignedPermutations(domainSize);
            end
        end


        %% Domain methods

        function b = eqv(self, x, y)
            b = isequal(x, y);
        end

        function h = hash(self, x)
            h = replab.Domain.hashVector(x + self.domainSize);
        end

        %% Monoid methods

        function z = compose(self, x, y)
            z = x(abs(y)).*sign(y);
        end

        %% Group methods

        function y = inverse(self, x)
            n = self.domainSize;
            y = zeros(1, n);
            xAbs = abs(x);
            y(xAbs) = 1:n;
            invFlip = xAbs(x < 0);
            y(invFlip) = -y(invFlip);
        end

        %% NiceFiniteGroup methods

        function res = hasSameParentAs(self, rhs)
            res = isa(rhs, 'replab.SignedPermutationGroup') && (self.parent.domainSize == rhs.parent.domainSize);
        end

        function p1 = niceMonomorphismImage(self, p)
            p1 = replab.SignedPermutations.toPermutation(p);
        end

        function grp = subgroup(self, generators, order)
        % Constructs a permutation subgroup from its generators
        %
        % Args:
        %   generators (row cell array): List of generators given as a permutations in a row cell array
        %   order (vpi, optional): Argument specifying the group order, if given can speed up computations
        %
        % Returns:
        %   +replab.SignedPermutationGroup: The constructed signed permutation subgroup
            if nargin < 3
                order = [];
            end
            grp = replab.SignedPermutationGroup(self.domainSize, generators, order, self.parent);
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
            G = replab.Permutations(self.domainSize).subgroup(newGenerators);
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

end
