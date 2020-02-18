classdef Rep < replab.Str
% Describes a finite dimensional representation of a compact group
%
% This class has mutable properties that correspond to information that can be computed and cached
% after the `+replab.Rep` instance is constructed, for example `.isUnitary` or `.isTrivial`.
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
        isUnitary % ({true, false, []}): Whether this representation is unitary
        isTrivial % ({true, false, []}): Whether this representation is trivial (not necessarily of dimension 1)
        isIrreducible % (true, false, []): Whether this representation is irreducible
        frobeniusSchurIndicator % (double or []): Value of the Frobenius-Schur indicator
        isDivisionAlgebraCanonical % ({true, false, []}): If the representation is real and irreducible, describes if its division algebra has canonical form
    end

    properties (SetAccess = protected)
        group     % (`+replab.CompactGroup`): Group being represented
        field     % ({'R', 'C'}): Vector space defined on real (R) or complex (C) field
        dimension % (integer): Representation dimension
    end

    properties (Access = protected)
        commutant_
        hermitianInvariant_
        decomposition_
    end

    methods % Abstract methods

        function rho = image_internal(self, g)
        % Returns the image of a group element, dense or sparse
        %
        % Args:
        %   g (element of `group`): Element being represented
        %
        % Returns:
        %   double(*,*): Image of the given element for this representation, may be sparse
            error('Abstract');
        end

    end

    methods

        %% Own methods

        function rho = image(self, g)
        % Returns the image of a group element
        %
        % Args:
        %   g (element of `group`): Element being represented
        %
        % Returns:
        %   double(*,*): Image of the given element for this representation
            rho = full(self.image_internal(g));
        end

        function rho = inverseImage_internal(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `group`): Element of which the inverse is represented
        %
        % Returns:
        %   double(*,*): Image of the inverse of the given element for this representation, may be sparse
            gInv = self.group.inverse(g);
            rho = self.image_internal(gInv);
        end

        function rho = inverseImage(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `group`): Element of which the inverse is represented
        %
        % Returns:
        %   double(*,*): Image of the inverse of the given element for this representation
            rho = full(self.inverseImage_internal(g));
        end

        %% Sampling

        function [rho rhoInverse] = sample(self)
        % Samples an element from the group and returns its image under the representation
        %
        % Optionally, the inverse of the image can be returned as well.
        %
        % Returns
        % -------
        %   rho: double(*,*)
        %     Image of the random group element
        %   rhoInverse: double(*,*)
        %     Inverse of image ``rho``
            g = self.group.sample;
            rho = self.image(g);
            if nargout > 1
                rhoInverse = self.image(self.group.inverse(g));
            end
        end

        %% Properties

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

        %% Derived vector spaces/algebras

        function e = equivariant(self, repC)
        % Returns the space of equivariant linear maps from this rep to another rep
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``self.image(g) * X = X * repC.image(g)``
        %
        %
        % Args:
        %   repC (replab.Rep): Representation on the source/column space
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            e = replab.makeEquivariant(self, repC, []);
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
            if isempty(self.commutant_)
                self.commutant_ = replab.makeEquivariant(self, self, 'commutant');
            end
            c = self.commutant_;
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
            if isempty(self.hermitianInvariant_)
                self.hermitianInvariant_ = replab.makeEquivariant(self, self.dual.conj, 'hermitian');
            end
            h = self.hermitianInvariant_;
        end

        %% Irreducible decomposition

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
            assert(isequal(self.isUnitary, true), 'Representation must be unitary');
            if isempty(self.decomposition_)
                dec = replab.irreducible.decomposition(self);
                self.decomposition_ = dec;%.nice;
            end
            I = self.decomposition_;
        end

        %% Str methods

        function s = headerStr(self)
            if self.isUnitary
                switch self.field
                  case 'R'
                    f = 'Orthogonal real';
                  case 'C'
                    f = 'Unitary complex';
                end
            else
                switch self.field
                  case 'R'
                    f = 'Real';
                  case 'C'
                    f = 'Complex';
                end
            end
            s = sprintf('%s representation of dimension %d', f, self.dimension);
        end
        %% Derived actions

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
        %   M (double(*,*)): Matrix acted upon
        %
        % Returns:
        %   double(*,*): The matrix ``self.image(g) * M``
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
        %   M (double(*,*)): Matrix acted upon
        %
        % Returns:
        %   double(*,*): The matrix ``M * self.inverseImage(g)``
            M = full(M * self.inverseImage_internal(g));
        end

        %% Derived representations

        function complexRep = complexification(self)
        % Returns the complexification of a real representation
        %
        % Returns:
        %   `+replab.Rep`: The complexification of this representation
        %
        % Raises:
        %   An error if this representation is already complex.
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.simplify(replab.rep.ComplexifiedRep(self));
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
            rep = replab.rep.simplify(replab.rep.DerivedRep(self, true, false, false));
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
            rep = replab.rep.simplify(replab.rep.DerivedRep(self, false, true, true));
        end

        function rep = blkdiag(varargin)
        % Direct sum of representations
        %
        % See `+replab.Rep.directSum`
            rep = replab.Rep.directSum(varargin);
        end

        function rep = kron(varargin)
        % Tensor product of representations
        %
        % See `+replab.Rep.tensor`
            rep = replab.Rep.tensor(varargin);
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
            rep = replab.Rep.tensor(reps);
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
            rep = replab.Rep.directSum(reps);
        end

        %% Manipulation of representation space

        function [A Ainv] = unitaryChangeOfBasis(self)
        % Returns the change of basis to a unitary representation
        %
        % Returns ``A`` and ``Ainv`` so that ``A * self.image(g) * Ainv`` is unitary.
        %
        % Returns
        % -------
        %   A: double(*,*)
        %     Change of basis matrix
        %   Ainv: double(*,*)
        %     Inverse of change of basis matrix
            if isequal(self.isUnitary, true)
                A = eye(self.dimension);
                Ainv = eye(self.dimension);
            else
                E = self.dual.conj.equivariant(self);
                X = E.project(eye(self.dimension));
                for i = 1:self.dimension
                    X(i,i) = real(X(i,i));
                end
                A = chol(X);
                Ainv = inv(A);
            end
        end

        function [newRep A Ainv] = unitarize(self)
        % Returns a unitary representation equivalent to this representation
        %
        % We have ``newRep.image(g) = A * self.image(g) * Ainv``.
        %
        %
        % Example:
        %   >>> S3 = replab.Permutations(3);
        %   >>> defRep = S3.definingRep.complexification;
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
        % Returns
        % -------
        %   newRep: `+replab.Rep`
        %     Unitary representation
        %   A: double(*,*)
        %     Change of basis matrix
        %   Ainv: double(*,*)
        %     Inverse of change of basis matrix
            [A Ainv] = self.unitaryChangeOfBasis;
            if isequal(self.isUnitary, true)
                newRep = self;
            else
                newRep = self.similar(A, Ainv);
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
        %   basis (double(dParent,dChild)): Basis of the subrepresentation
        %   embedding (double(dChild,dParent), optional): Map from the parent space to the subrepresentation
        % Returns:
        %   `+replab.SubRep`: Subrepresentation
            if nargin < 3
                dSub = size(basis, 2);
                rest = null(basis.');
                X = [basis rest];
                Xinv = inv(X);
                P = basis * Xinv(1:dSub, :);
                P1 = self.commutant.project(P);
                embedding = basis \ P1;
            end
            sub = replab.SubRep(self, basis, embedding, []);
        end

        function rep1 = similar(self, A, Ainv)
        % Returns a similar representation under a change of basis
        %
        % It returns a representation ``rep1`` such that
        %
        % ``rep1.image(g) = A * self.image(g) * Ainv``
        %
        % Args:
        %   A (double(*,*)): Change of basis matrix
        %   Ainv (double(*,*)): Inverse of the change of basis matrix
        %
        % Returns:
        %   `+replab.SimilarRep`: The similar representation
            rep1 = replab.SimilarRep(self, A, Ainv);
        end

    end

    methods (Static)

        function rep = directSum(reps)
        % Computes the direct sum of representations
        %
        % Args:
        %   reps (cell(1,*) of `+replab.Rep`): Representation of the same group over the same field
        %
        % Returns:
        %   `+replab.Rep`: Direct sum of the representations
            rep = replab.rep.DirectSumRep(reps);
        end

        function rep = tensor(reps)
        % Computes the tensor product of representations
        %
        % Args:
        %   reps (cell(1,*) of `+replab.Rep`): Representation of the same group over the same field
        %
        % Returns:
        %   `+replab.Rep`: Tensor product of the representations
            rep = replab.rep.TensorRep(reps);
        end

        function rep = lambda(group, field, dimension, isUnitary, irrepInfo, imageFun, inverseImageFun)
        % Creates a non unitary representation from an image function
        %
        % Args:
        %   group (replab.Group): Group represented
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   isUnitary (logical or []): Whether the representation is unitary, or ``[]`` if unknown
        %   irrepInfo (`+replab.+irreducible.Info` or []): Info about the irreducibility of this rep., or ``[]`` if unknown
        %   imageFun (function_handle): Function handle that returns an image matrix given a group element
        %   inverseImageFun (function_handle): Function handle that returns the inverse of the image
        %                                      matrix given a group element
        %
        % Returns:
        %   `+replab.Rep`: The constructed representation
            rep = replab.lambda.Rep(group, field, dimension, isUnitary, irrepInfo, imageFun, inverseImageFun);
        end

    end

end
