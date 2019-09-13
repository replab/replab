classdef Rep < replab.Str
% Describes a finite dimensional representation of a compact group
% 
% Only "image" needs to be implemented in principle.
%
% For optimization purposes, actions can also be specialized.
    
    properties (SetAccess = protected)
        group     % replab.CompactGroup: Group being represented
        field     % {'R', 'C'}: Vector space defined on real (R) or complex (C) field
        dimension % integer: Representation dimension
        isUnitary % {true, false, []}: Whether the representation is unitary
    end
    
    properties (Access = protected)
        commutant_ = [];
        decomposition_ = [];
    end
        
    methods % Abstract methods
        
        function rho = image(self, g)
        % Returns the image of a group element
        %
        % Args:
        %   g (element of `group`): Element being represented 
        %
        % Returns: 
        %   double matrix: Image of the given element for this representation
            error('Abstract');
        end
        
    end
    
    methods
        
        %% Own methods
        
        function rho = inverseImage(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `group`): Element of which the inverse is represented
        %
        % Returns:
        %   double matrix: Image of the inverse of the given element for this representation
            gInv = self.group.inverse(g);
            rho = self.image(gInv);
        end

        %% Sampling
        
        function [rho rhoInverse] = sample(self)
        % Samples an element from the group and returns its image under the representation
        %
        % Optionally, the inverse of the image can be returned as well.
        %
        % Returns
        % -------
        %   rho: double matrix
        %     Image of the random group element
        %   rhoInverse: double matrix
        %     Inverse of image `rho`
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
        % self.image(g) * X = X * repC.image(g)
        %
        %
        % Args:
        %   repC (replab.Rep): Representation on the source/column space
        %
        % Returns:
        %   replab.Equivariant: The equivariant vector space
            e = replab.EquivariantDispatch.instance.call(self, repC);
        end

        function c = commutant(self)
        % Returns the commutant of this representation
        %
        % This is the algebra of matrices that commute with the representation,
        % i.e. the vector space isomorphism to the equivariant space from this rep to this rep.
        %
        % For any g in G, we have ``rho(g) * X = X * rho(g)``.
        %
        % The computation is cached.
        %
        % Returns:
        %   replab.Commutant: The commutant algebra
            if isempty(self.commutant_)
                self.commutant_ = replab.Commutant(self);
            end
            c = self.commutant_;
        end
        
        %% Irreducible decomposition
        
        function I = decomposition(self)
        % Returns the irreducible decomposition of this representation
        %
        % Requires this representation to be unitary
        %
        % Returns:
        %   replab.Irreducible: The irreducible decomposition
        %
        % Raises:
        %   An error is this representation is not unitary.
            assert(~isempty(self.isUnitary) && self.isUnitary, 'Representation must be unitary');
            if isempty(self.decomposition_)
                dec = replab.rep.decomposition(self);
                self.decomposition_ = dec.nice;
            end
            I = self.decomposition_;
        end

        %% Str methods

        function s = headerStr(self)
            f = replab.str.field(self.field, 'Orthogonal real', 'Unitary complex');
            s = sprintf('%s representation of dimension %d', f, self.dimension);
        end
        
        %% Derived actions
        
        function M = matrixRowAction(self, g, M)
        % Computes the representation-matrix product
        %
        % This is a left action:
        %
        % self.matrixRowAction(g, self.matrixRowAction(h, M))
        %
        % is the same as
        %
        % self.matrixRowAction(self.group.compose(g, h), M)
        %
        % Args:
        %   g (`self.group` element): Group element acting
        %   M (double matrix): Matrix acted upon
        % 
        % Returns:
        %   double matrix: The matrix self.image(g) * M
            M = self.image(g) * M;
        end
        
        function M = matrixColAction(self, g, M)
        % Computes the matrix-representation product
        %
        % We multiply by the inverse of the image, so this stays a left action.
        %
        % self.matrixColAction(g, self.matrixColAction(h, M))
        %
        % is the same as
        %
        % self.matrixColAction(self.group.compose(g, h), M)

        % Args:
        %   g (`self.group` element): Group element acting
        %   M (double matrix): Matrix acted upon
        %
        % Returns:
        %   double matrix: The matrix M * self.inverseImage(g)
            M = M * self.inverseImage(g);
        end
        
        %% Derived representations
        
        function complexRep = complexification(self)
        % Returns the complexification of a real representation
        %
        % Raises:
        %   An error if this representation is already complex.
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.Rep.lambda(self.group, 'C', self.dimension, self.isUnitary, @(g) self.image(g), @(g) self.inverseImage(g));
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
        %   replab.Rep: The complex conjugate of this representation
            if self.overR
                rep = self;
            else
                rep = replab.Rep.lambda(self.group, 'C', self.dimension, self.isUnitary, @(g) conj(self.image(g)), @(g) conj(self.inverseImage(g)));
            end
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
            if self.isUnitary
                % If this is unitary, less drama to simply return the conjugate
                rep = self.conj;
            else
                imageFun = @(g) self.inverseImage(g).';
                inverseImageFun = @(g) self.image(g).';
                rep = replab.Rep.lambda(self.group, self.field, self.dimension, self.isUnitary, imageFun, inverseImageFun);
            end
        end
        
        function rep = blkdiag(varargin)
            rep = replab.Rep.directSum(varargin);            
        end
        
        function rep = kron(varargin)
            rep = replab.Rep.tensor(varargin);
        end
        
        %% Manipulation of representation space

        function sub = subRep(self, F, G)
        % Returns a subrepresentation of this representation
        %
        % Let V be the vector space of this representation, 
        % and W an invariant subspace of V. Let ``F: V -> W``
        % and ``G: W -> V`` be maps such that ``F G = id``,
        % and ``G rho(g) F`` is the subrepresentation.
        % 
        % Args:
        %   F (double matrix): map V -> W
        %   G (double matrix): map W -> V
        % Returns:
        %   replab.NUSubRep: Subrepresentation
            sub = replab.SubRepNU(self, F, G);
        end

        function sub = subRepUnitary(self, U)
        % Returns a unitary subrepresentation of this unitary representation
        %
        % It is described by the row basis vectors in U, such that
        % sub.image(g) = U * self.image(g) * U'
        %
        % U has dimension dim(subRepresentation) x self.dimension
        %
        % U needs to be orthogonal; if U is not orthonormal, the
        % basis vectors will be implicitly normalized
        %
        % Args:
        %   U: Orthogonal basis vectors stored as row vectors
        %
        % Returns:
        %   replab.SubRep: The subrepresentation
            assert(~isempty(self.isUnitary) && self.isUnitary, 'Representation must be unitary');
            sub = replab.SubRep(self, U);
        end

        function rep1 = leftConjugateUnitary(self, A)
        % Returns the (left) conjugation of this representation
        %
        % A must be a unitary matrix, and this representation must be unitary.
        %
        % It returns a representation rep1 such that
        %
        % rep1.image(g) = A * self.image(g) * inv(A)
        %
        % Args:
        %   A (double matrix): Conjugation matrix
        %
        % Returns:
        %   replab.ConjugateRep: The conjugated representation
            assert(~isempty(self.isUnitary) && self.isUnitary, 'Representation must be unitary');
            rep1 = replab.ConjugateRep(A, self);
        end

    end

    methods (Static)

        function rep = directSum(reps)
        % Computes the direct sum of representations
        %
        % Args:
        %   reps (row cell array of replab.Rep): Representation of the same group over the same field
        %
        % Returns:
        %   replab.Rep: Direct sum of the representations
            rep = replab.rep.DirectSumRep(reps);
        end

        function rep = tensor(reps)
        % Computes the tensor product of representations
        %
        % Args:
        %   reps (row cell array of replab.Rep): Representation of the same group over the same field
        %
        % Returns:
        %   replab.Rep: Tensor product of the representations
            rep = replab.rep.TensorRep(reps);
        end

        function rep = lambda(group, field, dimension, isUnitary, imageFun, inverseImageFun)
        % Creates a non unitary representation from an image function
        %
        % Args:
        %   group (replab.Group): Group represented
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %   imageFun (function_handle): Function handle that returns an image matrix given a group element
        %   inverseImageFun (function_handle, optional): Function handle that returns the inverse of the image
        %                                                matrix given a group element
        %
        % Returns:
        %   replab.Rep: The constructed representation
            if nargin < 6
                inverseImageFun = [];
            end
            rep = replab.lambda.Rep(group, field, dimension, isUnitary, imageFun, inverseImageFun);
        end

    end

end
