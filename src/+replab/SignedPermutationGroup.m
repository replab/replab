classdef SignedPermutationGroup < replab.gen.FiniteGroup
% A base class for all signed permutation groups

    properties (SetAccess = protected)
        domainSize % (integer): The integer $d$ when this group acts on $\{-d, .., -1, 1, .., d\}$
    end

    methods

        function self = SignedPermutationGroup(domainSize, generators, varargin)
        % Constructs a signed permutation group
        %
        % Additional keyword arguments (``type``, ``nice`` and ``niceIsomorphism``) are used internally
        % but are not part of the public API.
        %
        % Args:
        %   domainSize (integer): Domain size of this permutation group
        %   generators (cell(1,\*) of integer(1,\*)): Group generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring): Names of the generators
        %   order (vpi, optional): Order of the group
        %   relators (cell(1,\*) of charstring): Relators given either in word or letter format
            args = struct('type', []);
            [args, rest] = replab.util.populateStruct(args, varargin);
            if isempty(args.type)
                type = replab.signed.FiniteGroupType.make(domainSize);
            else
                type = args.type;
            end
            self@replab.gen.FiniteGroup(type, generators, rest{:});
            self.domainSize = domainSize;
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.SignedPermutationGroupLaws(self);
        end

        % Domain

        function b = eqv(self, x, y)
            b = all(x == y);
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

    end

    methods % Signed permutation methods

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

    end

    methods % Actions

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

    end

    methods % Representation construction

        function rho = naturalRep(self)
        % Natural representation on R^d of signed permutations on integers -d..-1, 1..d
            rho = self.signedPermutationRep(self.domainSize, 'preimages', self.generators, 'images', self.generators);
        end

    end

    methods (Static)

        function G = of(varargin)
        % Constructs a nontrivial signed permutation group from the given generators
        %
        % Example:
        %   >>> s = [-1 -2 -3 -4];
        %   >>> i = [2 -1 4 -3];
        %   >>> j = [3 -4 -1 2];
        %   >>> Q = replab.SignedPermutationGroup.of(s, i, j);
        %   >>> Q.order == 8
        %     1
        %
        % The generators of the group can be named by preceding them all by a charstring:
        %
        % Example:
        %   >>> s = [-1 -2 -3 -4];
        %   >>> i = [2 -1 4 -3];
        %   >>> j = [3 -4 -1 2];
        %   >>> Q = replab.SignedPermutationGroup.of('s', s, 'i', i, 'j', j);
        %   >>> Q.order == 8
        %     1
        %
        % This method cannot construct trivial groups without any generators.
        % In that case, use the constructor:
        %
        % Example:
        %   >>> generators = {};
        %   >>> domainSize = 4;
        %   >>> G = replab.SignedPermutationGroup(domainSize, generators);
        %   >>> G.order
        %       1
        %
        % Args:
        %   varargin (cell(1,\*) of signed permutations): Group generators, pos
        %
        % Returns:
        %   `.SignedPermutationGroup`: The signed permutation group given as the closure of the generators
            assert(nargin > 0, 'Must be called with at least one generator');
            mask = cellfun(@ischar, varargin);
            if mask(1) % named generators
                assert(mod(nargin, 2) == 0);
                assert(all(mask(1:2:end)));
                assert(all(~mask(2:2:end)));
                names = varargin(1:2:end);
                generators = varargin(2:2:end);
                domainSize = length(generators{1});
                G = replab.SignedPermutationGroup(domainSize, generators, 'generatorNames', names);
            else
                generators = varargin;
                domainSize = length(generators{1});
                G = replab.SignedPermutationGroup(domainSize, generators);
            end
        end

        function G = signedSymmetric(n)
        % Creates the signed symmetric group of a given degree
        %
        % Args:
        %   n (integer): Degree of the signed symmetric group
        %
        % Returns:
        %   `.SignedPermutationGroup`: The signed symmetric group of degree ``n``
            T = replab.signed.FiniteGroupType.make(n);
            G = T.isomorphism.source;
        end

    end

end
