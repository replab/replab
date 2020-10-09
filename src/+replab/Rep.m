classdef Rep < replab.Obj
% Describes a finite dimensional representation of a compact group
%
% This class has mutable properties that correspond to information that can be computed and cached
% after the `+replab.Rep` instance is constructed, for example `.isUnitary` or `.trivialDimension`.
%
% By convention, those properties are not passed to constructors of subclasses of `+replab.Rep`,
% but are instead their default value (``[]`` signifying unknown) is overwritten after construction.
%
% Notes:
%   While we do not expect users to implement their own subclass of `~+replab.Rep`, to do
%   so only the map from group elements to matrix images need to be implemented.
%
%   Internally RepLAB can use sparse matrices; however, we take care to not return sparse
%   matrices in the user API, as doing so can lead to spurious breakages in user code
%   (some MATLAB/Octave functions do not work with sparse arguments).
%
%   Thus, internally, the methods `.image_internal` and `.inverseImage_internal` are used,
%   and the user facing methods `.image` and `.inverseImage` simply process the output of
%   the internal methods by converting them, if necessary, to dense matrices
%
%   This class implements also row and column matrix actions, methods which can be overloaded
%   for performance.

    properties
        % see `+replab.+rep.copyProperties`
        isUnitary % ({true, false, []}): Whether this representation is unitary
        trivialDimension % (integer or []): Dimension of the trivial subrepresentation in this representation
        isIrreducible % (true, false, []): Whether this representation is irreducible
        frobeniusSchurIndicator % (integer or []): Value of the non-approximated Frobenius-Schur indicator
        isDivisionAlgebraCanonical % ({true, false, []}): If the representation is real and irreducible and the Frobenius-Schur indicator is not 1, means that the images encode the complex or quaternion division algebras in the RepLAB canonical form
    end

    properties
        nearEps2_ % (double): Upper bound on ``||rep(g) rep(h) - rep(gh)||2 <= nearEps2``
        nearEpsF_ % (double): Upper bound on ``||rep(g) rep(h) - rep(gh)||F <= nearEpsF``
        approxEps2_ % (double): There exists an exact representation ``repE`` such that ``||rep(g) - repE(g)||2 <= approxEps2``
        approxEpsF_ % (double): There exists an exact representation ``repE`` such that ``||rep(g) - repE(g)||2 <= approxEpsF``
    end

    methods

        function computeNearEps(self)
            if isa(self.group, 'replab.FiniteGroup')
                els = self.group.elements.toCell;
                n = length(els);
                for i = 1:n
                    for j = 1:n
                        g = els{i};
                        h = els{j};
                        gh = self.group.compose(g, h);
                        repgh1 = self.image(g) * self.image(h);
                        repgh2 = self.image(gh);
                        diff = repgh1 - repgh2;
                        self.nearEps2_ = norm(diff, 2);
                        self.nearEpsF_ = norm(diff, 'fro');
                    end
                end
            end
        end

        function e = nearEps2(self)
            if isempty(self.nearEps2_)

                e = NaN;
            else
                e = self.nearEps2_;
            end
        end

        function e = nearEpsF(self)
            if isempty(self.nearEpsF_)
                e = NaN;
            else
                e = self.nearEpsF_;
            end
        end

        function e = approxEps2(self)
            if isempty(self.approxEps2_)
                e = NaN;
            else
                e = self.approxEps2_;
            end
        end

        function e = approxEpsF(self)
            if isempty(self.approxEpsF_)
                e = NaN;
            else
                e = self.approxEpsF_;
            end
        end

    end

    properties (SetAccess = protected)
        group     % (`+replab.CompactGroup`): Group being represented
        field     % ({'R', 'C'}): Vector space defined on real (R) or complex (C) field
        dimension % (integer): Representation dimension
    end

    methods % Representation methods

        % Abstract method

        function [rho error2 errorF] = image_error(self, g)
        % Returns the image of a group element, dense or sparse, along with an error estimate
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %
        % Returns
        % -------
        %   rho:
        %     double(\*,\*): Image, may be sparse
        %   error2:
        %     double: Error estimation (2-norm)
        %   errorF:
        %     double: Error estimation (Frobenius norm)
        end

        function rho = image_internal(self, g)
        % Returns the image of a group element, dense or sparse
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %
        % Returns:
        %   double(\*,\*): Image of the given element for this representation, may be sparse
            error('Abstract');
        end

        % Methods with default implementations

        function rho = image(self, g)
        % Returns the image of a group element
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %
        % Returns:
        %   double(\*,\*): Image of the given element for this representation
            rho = full(self.image_internal(g));
        end

        function rho = inverseImage_internal(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `.group`): Element of which the inverse is represented
        %
        % Returns:
        %   double(\*,\*): Image of the inverse of the given element for this representation, may be sparse
            gInv = self.group.inverse(g);
            rho = self.image_internal(gInv);
        end

        function rho = inverseImage(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `.group`): Element of which the inverse is represented
        %
        % Returns:
        %   double(\*,\*): Image of the inverse of the given element for this representation
            rho = full(self.inverseImage_internal(g));
        end

        function [rho rhoInverse] = sample(self)
        % Samples an element from the group and returns its image under the representation
        %
        % Optionally, the inverse of the image can be returned as well.
        %
        % Returns
        % -------
        %   rho: double(\*,\*)
        %     Image of the random group element
        %   rhoInverse: double(\*,\*)
        %     Inverse of image ``rho``
            g = self.group.sample;
            rho = self.image(g);
            if nargout > 1
                rhoInverse = self.image(self.group.inverse(g));
            end
        end

        % Properties

        function b = overR(self)
        % Returns whether this representation is defined over the real field
        %
        % Returns:
        %   logical: True if this representation is defined over the real field
            b = isequal(self.field, 'R');
        end

        function b = overC(self)
        % Returns whether this representation is defined over the complex field
        %
        % Returns:
        %   logical: True if this representation is defined over the complex field
            b = isequal(self.field, 'C');
        end

    end

    methods % Computed properties

        function K = kernel(self)
        % Returns the kernel of the given representation
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> K = S3.signRep.kernel;
        %   >>> K == replab.PermutationGroup.alternating(3)
        %       1
        %
        % Returns:
        %   `+replab.FiniteGroup`: The group ``K`` such that ``rho.image(k) == id`` for all ``k`` in ``K``
            K = self.cached('kernel', @() self.computeKernel);
        end

        function K = computeKernel(self)
            assert(isa(self.group, 'replab.FiniteGroup'));
            % TODO error estimation: take in account the uncertainty on computed images
            C = self.group.conjugacyClasses.classes;
            % for a character, we have chi(g) == chi(id) only if rho(g) == eye(d)
            % what is the maximal value of real(chi(g)) for chi(g) ~= chi(id)?
            % write chi(g) = sum(lambda(g)) where lambda(g) = eig(chi(g))
            % the worst case is when lambda(g)(2:d) = 1, and lambda(g)(1) is something else
            % now lambda(g)(1) ~= 1; we know that lambda(g)(1) is a root of unity, and that
            % lambda(g)(1)^group.exponent == 1
            % thus, the maximum real part is when lambda(g)(1) == exp(2*pi/group.exponent)
            % which has real part cos(2*pi/group.exponent)
            maxNTChi = self.dimension-cos(2*pi/double(self.group.exponent)); % maximum value for nontrivial character
            chi = cellfun(@(c) trace(self.image_internal(c.representative)), C);
            mask = real(chi) > maxNTChi;
            K = self.group.subgroup(cellfun(@(c) c.representative, C(mask), 'uniform', 0));
            K = self.group.normalClosure(K);
        end

    end

    methods % Derived vector spaces/algebras

        function e = equivariantTo(self, repR)
        % Returns the space of equivariant linear maps from this rep to another rep
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * self.image(g) = repR.image(g) * X``
        %
        %
        % Args:
        %   repR (`+replab.Rep`): Representation on the target/row space
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            e = replab.Equivariant.make(self, repR, '');
        end

        function c = computeCommutant(self)
            c = replab.Equivariant.make(self, self, 'commutant');
        end

        function c = commutant(self)
        % Returns the commutant of this representation
        %
        % This is the algebra of matrices that commute with the representation,
        % i.e. the vector space isomorphism to the equivariant space from this rep to this rep.
        %
        % For any ``g in G``, we have ``rho(g) * X = X * rho(g)``.
        %
        % The computation is cached.
        %
        % Returns:
        %   `+replab.Equivariant`: The commutant algebra represented as an equivariant space
            c = self.cached('commutant', @() self.computeCommutant);
        end

        function h = computeHermitianInvariant(self)
            h = replab.Equivariant.make(self.dual.conj, self, 'hermitian');
        end

        function h = hermitianInvariant(self)
        % Returns the Hermitian invariant space of this representation
        %
        % This is the space of Hermitian matrices that are invariant under this representation
        % i.e.
        %
        % for any g in G, we have ``X = rho(g) * X * rho(g^-1)'``
        %
        % Returns:
        %   `+replab.Equivariant`: The space of Hermitian invariant matrices
            h = self.cached('hermitianInvariant', @() self.computeHermitianInvariant);
        end

        function t = computeTrivialSpace(self)
            tRep = self.group.trivialRep(self.field, self.dimension);
            t = replab.Equivariant.make(tRep, self, 'trivial');
        end

        function t = trivialSpace(self)
        % Returns an equivariant space from a trivial representation to this representation
        %
        % The trivial representation has the same dimension as this representation
            t = self.cached('trivialSpace', @() self.computeTrivialSpace);
        end

    end

    methods % Irreducible decomposition

        function I = decomposition(self)
        % Returns the irreducible decomposition of this representation
        %
        % Requires this representation to be unitary
        %
        % Returns:
        %   `+replab.Irreducible`: The irreducible decomposition
        %
        % Raises:
        %   An error is this representation is not unitary.
            I = self.cached('decomposition', @() self.computeDecomposition);
        end

        function dec = computeDecomposition(self)
            dec = replab.irreducible.decomposition(self);
            if dec.nComponents == 1 && dec.components{1}.multiplicity == 1
                assert(~isequal(self.isIrreducible, false));
                self.isIrreducible = true;
                if isequal(dec.basis, eye(self.dimension))
                    replab.rep.copyProperties(dec, self);
                end
            end
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            p = {};
            if isequal(self.isUnitary, true)
                if self.overR
                    p{1,end+1} = 'orthogonal';
                else
                    p{1,end+1} = 'unitary';
                end
            elseif isequal(self.isUnitary, false)
                if self.overR
                    p{1,end+1} = 'nonorthogonal';
                else
                    p{1,end+1} = 'nonunitary';
                end
            else
                if self.overR
                    p{1,end+1} = 'real';
                else
                    p{1,end+1} = 'complex';
                end
            end
            if isequal(self.trivialDimension, self.dimension)
                p{1,end+1} = 'trivial';
            elseif isequal(self.trivialDimension, 0)
                p{1,end+1} = 'fully nontrivial';
            elseif ~isempty(self.trivialDimension)
                p{1,end+1} = 'nontrivial';
            end
            if isequal(self.isIrreducible, true)
                p{1,end+1} = 'irreducible';
            elseif isequal(self.isIrreducible, true)
                p{1,end+1} = 'reducible';
            end
            if ~isequal(self.frobeniusSchurIndicator, []) && self.overR
                switch self.frobeniusSchurIndicator
                  case 1
                    p{1,end+1} = 'real-type';
                  case 0
                    p{1,end+1} = 'complex-type';
                  case -2
                    p{1,end+1} = 'quaternion-type';
                  otherwise
                    % do nothing
                end
            end
            switch class(self)
              case 'replab.SubRep'
                p{1,end+1} = 'subrepresentation';
              case 'replab.rep.DerivedRep'
                p{1,end+1} = 'derived representation';
                if self.conjugate
                    p{1,end+1} = '(conjugate)';
                end
                if self.inverse
                    p{1,end+1} = '(inverse)';
                end
                if self.transpose
                    p{1,end+1} = '(transpose)';
                end
              case 'replab.rep.DirectSumRep'
                p{1,end+1} = 'direct sum representation';
              case 'replab.rep.TensorRep'
                p{1,end+1} = 'tensor representation';
              case 'replab.rep.ComplexifiedRep'
                p{1,end+1} = 'complexified representation';
              case 'replab.RepByImages'
                p{1,end+1} = 'representation by images';
              otherwise
                p{1,end+1} = 'representation';
            end
            p{1} = replab.str.capitalize(p{1});
            s = strjoin(p, ' ');
        end

        % Obj

        function l = laws(self)
            l = replab.laws.RepLaws(self);
        end

    end

    methods % Derived actions

        function M = matrixRowAction(self, g, M)
        % Computes the representation-matrix product
        %
        % This is a left action:
        %
        % ``self.matrixRowAction(g, self.matrixRowAction(h, M))``
        %
        % is the same as
        %
        % ``self.matrixRowAction(self.group.compose(g, h), M)``
        %
        % Args:
        %   g (`group` element): Group element acting
        %   M (double(\*,\*)): Matrix acted upon
        %
        % Returns:
        %   double(\*,\*): The matrix ``self.image(g) * M``
            M = full(self.image_internal(g) * M);
        end

        function M = matrixColAction(self, g, M)
        % Computes the matrix-representation product
        %
        % We multiply by the inverse of the image, so this stays a left action.
        %
        % ``self.matrixColAction(g, self.matrixColAction(h, M))``
        %
        % is the same as
        %
        % ``self.matrixColAction(self.group.compose(g, h), M)``
        %
        % Args:
        %   g (`group` element): Group element acting
        %   M (double(\*,\*)): Matrix acted upon
        %
        % Returns:
        %   double(\*,\*): The matrix ``M * self.inverseImage(g)``
            M = full(M * self.inverseImage_internal(g));
        end

    end

    methods % Morphism composition

        function res = imap(self, f)
        % Maps the representation under an isomorphism
        %
        % Args:
        %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
        %
        % Returns:
        %   `.Rep`: Representation satisfying ``newRep.group.isSubgroup(f.target)``.
            if self.group.order < f.source.order
                f = f.restrictedSource(self.group);
            end
            res = replab.rep.CompositionRep(f.inverse, self);
            res.isUnitary = self.isUnitary;
            res.trivialDimension = self.trivialDimension;
            res.isIrreducible = self.isIrreducible;
            res.frobeniusSchurIndicator = self.frobeniusSchurIndicator;
            res.isDivisionAlgebraCanonical = self.isDivisionAlgebraCanonical;
        end

        function res = compose(self, applyFirst)
        % Composition of a representation with a morphism, with the morphism applied first
        %
        % Args:
        %   applyFirst (`.Morphism`): Morphism from a finite group
        %
        % Returns:
        %   `.Rep`: Representation
            res = replab.rep.CompositionRep(self, applyFirst);
        end

    end

    methods % Derived representations

        function rep1 = contramap(self, morphism)
        % Returns the representation composed with the given morphism applied first
        %
        % Args:
        %   morphism (`+replab.FiniteMorphism`): Morphism of finite groups such that ``morphism.target == self.group``
        %
        % Returns:
        %   `+replab.Rep`: Representation on the finite group ``morphism.source``
            assert(self.group == morphism.target);
            rep1 = replab.Rep.lambda(morphism.source, self.field, self.dimension, @(g) rep.image_internal(morphism.image(g)), @(g) replab.inverseImage_internal(morphism.image(g)));
        end

        function rep1 = restrictedTo(self, subgroup)
        % Returns the restricted representation to the given subgroup
        %
        % Args:
        %   subgroup (`.CompactGroup`): Subgroup to restrict the representation to, must be a subgroup of `.group`
            rep1 = replab.Rep.lambda(subgroup, self.field, self.dimension, @(g) self.image(g), @(g) self.inverseImage(g));
        end

        function rep1 = simplify(self)
        % Returns a representation identical to this, but possibly with its composition simplified
        %
        % Returns:
        %   `+replab.Rep`: A possibly simplified representation
            rep1 = replab.rep.simplify(self);
        end

        function complexRep = complexification(self)
        % Returns the complexification of a real representation
        %
        % Returns:
        %   `+replab.Rep`: The complexification of this representation
        %
        % Raises:
        %   An error if this representation is already complex.
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.ComplexifiedRep(self);
        end

        function rep = conj(self)
        % Returns the complex conjugate representation of this representation
        %
        % See https://en.wikipedia.org/wiki/Complex_conjugate_representation
        %
        % It obeys ``rep.conj.image(g) = conj(rep.image(g))``
        %
        % If this representation is real, it is returned unchanged.
        %
        % Returns:
        %   `+replab.Rep`: The complex conjugate of this representation
            rep = replab.rep.DerivedRep(self, true, false, false);
        end

        function rep = dual(self)
        % Returns the dual representation of this representation
        %
        % See https://en.wikipedia.org/wiki/Dual_representation
        %
        % It obeys ``rep.dual.image(g) = rep.inverseImage(g).'``
        %
        % Returns:
        %   replab.Rep: The dual representation
            rep = replab.rep.DerivedRep(self, false, true, true);
        end

        function rep = blkdiag(varargin)
        % Direct sum of representations
        %
        % See `+replab.CompactGroup.directSumRep`
            self = varargin{1};
            rep = self.group.directSumRep(self.field, varargin);
        end

        function rep = kron(varargin)
        % Tensor product of representations
        %
        % See `+replab.CompactGroup.tensorRep`
            self = varargin{1};
            rep = self.group.tensorRep(self.field, varargin);
        end

        function rep = tensorPower(self, n)
        % Returns a tensor power of this representation
        %
        % Args:
        %   n (integer): Exponent of the tensor power
        %
        % Returns:
        %   `+replab.Rep`: The tensor power representation
            reps = arrayfun(@(i) self, 1:n, 'uniform', 0);
            rep = self.group.tensorRep(self.field, reps);
        end

        function rep = directSumOfCopies(self, n)
        % Returns a direct sum of copies of this representation
        %
        % Args:
        %   n (integer): Number of copies
        %
        % Returns:
        %   `+replab.Rep`: The direct sum representation
            reps = arrayfun(@(i) self, 1:n, 'uniform', 0);
            rep = self.group.directSum(self.field, reps);
        end

        %% Manipulation of representation space

        function res = computeUnitarize(self)
            if isequal(self.isUnitary, true)
                res = replab.SimilarRep.identical(self);
            else
                [A Ainv] = self.unitaryChangeOfBasis;
                res = self.similarRep(A, Ainv);
                res.isUnitary = true;
            end
        end

        function res = unitarize(self)
        % Returns a unitary representation equivalent to this representation
        %
        % The returned representation is of type `.SimilarRep`, from which
        % the change of basis matrix can be obtained.
        %
        % If the representation is already unitary, the returned `.SimilarRep`
        % has the identity matrix as a change of basis.
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> defRep = S3.naturalRep.complexification;
        %   >>> C = randn(3,3) + 1i * rand(3,3);
        %   >>> nonUnitaryRep = defRep.subRep(C, inv(C));
        %   >>> unitaryRep = nonUnitaryRep.unitarize;
        %   >>> U = unitaryRep.sample;
        %   >>> norm(U*U' - eye(3)) < 1e-10
        %      ans =
        %       logical
        %       1
        %
        %
        % Returns:
        %   `+replab.SimilarRep`: Unitary similar representation
            res = self.cached('unitarize', @() self.computeUnitarize);
        end

        function [sub1 sub2] = maschke(self, basis1, embedding1)
        % Given a basis of a subrepresentation, returns two complementary subrepresentations
        %
        % Note that we optimize special cases when the representation and/or the basis is
        % unitary
        %
        % Args:
        %   basis1 (double(dParent,dSub1)): Basis of the first subrepresentation
        %   embedding1 (double(dChild,dParent), optional): Map from the parent space to the subrepresentation
        %
        % Returns
        % -------
        %   sub1: `replab.SubRep`
        %     First subrepresentation
        %   sub2: `replab.SubRep`
        %     Second subrepresentation
            assert(size(basis1, 1) == self.dimension);
            d1 = size(basis1, 2);
            if nargin > 2
                assert(size(embedding1, 1) == d1);
                assert(size(embedding1, 2) == self.dimension);
            end
            rest = null(basis1.');
            if isequal(self.isUnitary, true)
                if nargin < 3
                    BB = basis1'*basis1;
                    BB = (BB+BB')/2;
                    embedding1 = BB \ basis1';
                end
                sub1 = replab.SubRep(self, basis1, embedding1);
                % for the second subrepresentation, we take the orthogonal complement from 'null'
                sub2 = replab.SubRep(self, rest, rest');
            else
                X = [basis1 rest];
                Xinv = inv(X);
                P = basis1 * Xinv(1:d1, :);
                P1 = self.commutant.project(P);
                embedding1 = basis1 \ P1;
                P2 = eye(self.dimension) - P1;
                basis2 = orth(P2);
                embedding2 = basis2 \ P2;
                sub1 = replab.SubRep(self, basis1, embedding1);
                sub2 = replab.SubRep(self, basis2, embedding2);
            end
        end

        function sub = subRep(self, basis, embedding)
        % Returns a subrepresentation of this representation
        %
        % The subrepresentation is defined by its basis in the parent representation; to compute
        % images, an embedding map can be provided.
        %
        % While the ``basis`` represents in essence a map from the subrepresentation to parent representation,
        % the embedding map is a map from the parent representation to the subrepresentation.
        %
        % The embedding map is not uniquely defined, for example when the subrepresentation contains irreducible
        % representations that have multiplicites outside the subrepresentation space.
        %
        % However, all variants of the embedding map provide identical results when computing images of the
        % subrepresentation.
        %
        % If the embedding is not provided, one is obtained by a trick based on Maschke theorem.
        % Args:
        %   basis (double(dParent,dChild), may be sparse): Basis of the subrepresentation
        %   embedding (double(dChild,dParent), may be sparse, optional): Map from the parent space to the subrepresentation
        % Returns:
        %   `+replab.SubRep`: Subrepresentation
            if nargin < 3
                if isequal(self.isUnitary, true)
                    % optimization for unitary parent representations
                    embedding = (basis'*basis) \ basis';
                    % there is an optimization that cannot be made yet:
                    % if basis'*basis is identity, then embedding = basis'
                    % works
                else
                    dSub = size(basis, 2);
                    rest = null(basis.');
                    X = [basis rest];
                    Xinv = inv(X);
                    P = basis * Xinv(1:dSub, :);
                    P1 = self.commutant.project(P);
                    embedding = basis \ P1;
                end
            end
            sub = replab.SubRep(self, basis, embedding);
        end

        function irreps = splitIntoIrreducibles(self, context)
        % Decomposes fully the given representation into subrepresentations
        %
        % Returns a list of irreducible representations, where trivial subrepresentations
        % have been identified
        %
        % If this representation is irreducible, it will set its `~+replab.Rep.isIrreducible`.
        %
        % Args:
        %   context (`+replab.Context`, optional): Sampling context to use
        %
        % Returns:
        %   cell(1,\*) of `+replab.SubRep`: irreducible subrepresentations
            if nargin < 2
                context = replab.Context.make;
            end
            d = self.dimension;
            start = replab.SubRep.fullSubRep(self);
            todo = {start};
            irreps = cell(1, 0);
            while ~isempty(todo)
                h = todo{1};
                if isequal(h.isIrreducible, true)
                    % head of list is irreducible, remove it
                    irreps{1,end+1} = h;
                    todo = todo(2:end);
                else
                    res = replab.irreducible.split(h, context);
                    res1 = cellfun(@(sub) replab.rep.collapse(sub), res, 'uniform', 0);
                    todo = horzcat(todo(2:end), res1);
                end
            end
            if nargin < 2
                context.close;
            end
        end

        function rep1 = similarRep(self, A, Ainv)
        % Returns a similar representation under a change of basis
        %
        % It returns a representation ``rep1`` such that
        %
        % ``rep1.image(g) = A * self.image(g) * Ainv``
        %
        % Args:
        %   A (double(\*,\*)): Change of basis matrix
        %   Ainv (double(\*,\*)): Inverse of the change of basis matrix
        %
        % Returns:
        %   `+replab.SimilarRep`: The similar representation
            rep1 = replab.SimilarRep(self, A, Ainv);
        end

    end

    methods (Access = protected)

        function [A Ainv] = unitaryChangeOfBasis(self)
        % Returns the change of basis to a unitary representation
        %
        % Returns ``A`` and ``Ainv`` so that ``A * self.image(g) * Ainv`` is unitary.
        %
        % Returns
        % -------
        %   A: double(\*,\*)
        %     Change of basis matrix
        %   Ainv: double(\*,\*)
        %     Inverse of change of basis matrix
            if isequal(self.isUnitary, true)
                A = eye(self.dimension);
                Ainv = eye(self.dimension);
            else
                X = self.hermitianInvariant.project(eye(self.dimension));
                for i = 1:self.dimension
                    X(i,i) = real(X(i,i));
                end
                Ainv = chol(X, 'lower');
                A = inv(Ainv);
            end
        end

    end

    methods (Static)

        function rep = lambda(group, field, dimension, image_internalFun, inverseImage_internalFun)
        % Creates a non unitary representation from an image function
        %
        % Args:
        %   group (replab.Group): Group represented
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   image_internalFun (function_handle): Function handle that returns an image matrix given a group element
        %   inverseImage_internalFun (function_handle): Function handle that returns the inverse of the image
        %                                               matrix given a group element
        %
        % Returns:
        %   `+replab.Rep`: The constructed representation
            rep = replab.lambda.Rep(group, field, dimension, image_internalFun, inverseImage_internalFun);
        end

    end

end
