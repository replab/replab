classdef Rep < replab.Obj
% Describes a finite dimensional representation of a compact group
%
% This class has mutable properties that correspond to information that can be computed and cached
% after the `+replab.Rep` instance is constructed, for example `.isUnitary` or `.trivialDimension`.
%
% Notes:
%   While we do not expect users to implement their own subclass of `~+replab.Rep`, to do
%   so only the map from group elements to matrix images need to be implemented.
%
%   Either the override the `.image_double` method if the user implementation provides floating-point
%   images, or both the `.isExact` and `.image_exact` methods.
%
%   This class implements also extra methods about the action of the representation on matrices,
%   methods which can be overloaded for performance.

    properties (SetAccess = protected)
        group     % (`+replab.CompactGroup`): Group being represented
        field     % ({'R', 'C'}): Vector space defined on real (R) or complex (C) field
        dimension % (integer): Representation dimension
        isUnitary % (logical): Whether the representation is unitary
    end

    methods

        function self = Rep(group, field, dimension, varargin)
        % Constructs a finite dimension representation
        %
        % Args:
        %   group (`+replab.CompactGroup`): Group being represented
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %
        % Keyword Args:
        %   isUnitary (logical): Whether the representation is unitary, mandatory argument
        %   isIrreducible (logical or ``[]``, optional): Whether this representation is irreducible, default ``[]``
        %   trivialDimension (integer or ``[]``, optional): Dimension of the trivial subrepresentation, default ``[]``
        %   frobeniusSchurIndicator % (integer or ``[]``, optional): Exact value of the Frobenius-Schur indicator, default ``[]``
        %   isDivisionAlgebraCanonical % (logical or ``[]``): If the representation is real and irreducible and the Frobenius-Schur indicator is not 1, means that the images encode the complex or quaternion division algebras in the RepLAB canonical form
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            args = struct('isUnitary', [], 'isIrreducible', [], 'trivialDimension', [], 'frobeniusSchurIndicator', [], 'isDivisionAlgebraCanonical', []);
            args = replab.util.populateStruct(args, varargin);
            assert(~isempty(args.isUnitary), 'The isUnitary keyword parameter must be provided.');
            self.isUnitary = args.isUnitary;
            if ~isempty(args.trivialDimension)
                self.cache('trivialDimension', args.trivialDimension, 'error');
            end
            if ~isempty(args.isIrreducible)
                self.cache('isIrreducible', args.isIrreducible, 'error');
            end
            if ~isempty(args.frobeniusSchurIndicator)
                self.cache('frobeniusSchurIndicator', args.frobeniusSchurIndicator, 'error');
            end
            if ~isempty(args.isDivisionAlgebraCanonical)
                self.cache('isDivisionAlgebraCanonical', args.isDivisionAlgebraCanonical, 'error');
            end
        end

    end

    methods (Access = protected) % Protected methods

        function e = computeErrorBound(self)
            error('Abstract');
        end

        function rho = image_double(self, g)
            rho = double(self.image_exact(g));
        end

        function rho = image_exact(self, g)
            error('Exact images not implemented');
        end

    end

    methods % Image computation

        function rho = image(self, g, type)
        % Returns the image of a group element
        %
        % Raises:
        %   An error if ``type`` is ``'exact'`` and ``self.isExact`` is false.
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %   type ({'double', 'exact'}, optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Image of the given element for this representation
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'double'
                rho = self.image_double(g);
              case 'exact'
                rho = self.image_exact(g);
              otherwise
                error('Type must be either double or exact');
            end
        end

        function b = isExact(self)
        % Returns whether this representation can compute exact images
        %
        % Returns:
        %   logical: True if the call ``self.image(g, 'exact')`` always succeeds
            b = false;
        end

        function rho = inverseImage(self, g, type)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `.group`): Element of which the inverse is represented
        %   type ({'double', 'exact'}, optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Image of the inverse of the given element for this representation
            if nargin < 3
                type = 'double';
            end
            if self.isUnitary
                rho = self.image(g, type);
                rho = rho';
            else
                rho = self.image(self.group.inverse(g), type);
            end
        end

        function [rho rhoInverse] = sample(self, type)
        % Samples an element from the group and returns its image under the representation
        %
        % Args:
        %   type ({'double', 'exact'}, optional): Type of the returned value, default: 'double'
        %
        % Optionally, the inverse of the image can be returned as well.
        %
        % Returns
        % -------
        %   rho: double(\*,\*)
        %     Image of the random group element
        %   rhoInverse: double(\*,\*)
        %     Inverse of image ``rho``
            if nargin < 2
                type = 'double';
            end
            g = self.group.sample;
            rho = self.image(g, type);
            if nargout > 1
                rhoInverse = self.inverseImage(g, type);
            end
        end

    end

    methods (Access = protected)

        function c = computeConditionNumberEstimate(self)
            if self.isUnitary
                c = 1;
            else
                simRep = self.unitarize;
                c = cond(simRep.A_internal);
            end
        end

        function b = inKernel(self, g)
        % Returns whether the given group element is in the kernel of the representation
        %
        % Args:
        %   g (group element): Group element to test
        %
        % Returns:
        %   logical: True if the group element is in the kernel
            if self.isExact
                chi = trace(self.image(g, 'exact'));
                b = chi == self.dimension;
            else
                % for a character, we have chi(g) == chi(id) only if rho(g) == eye(d)
                % what is the maximal value of real(chi(g)) for chi(g) ~= chi(id)?
                % write chi(g) = sum(lambda(g)) where lambda(g) = eig(chi(g))
                % the worst case is when lambda(g)(2:d) = 1, and lambda(g)(1) is something else
                % now lambda(g)(1) ~= 1; we know that lambda(g)(1) is a root of unity, and that
                % lambda(g)(1)^group.exponent == 1
                % thus, the maximum real part is when lambda(g)(1) == exp(2*pi/group.exponent)
                % which has real part cos(2*pi/group.exponent)
                nonTrivialUB = self.dimension-cos(2*pi/double(self.group.exponent)); % maximum for nontrivial character
                chi = real(trace(self.image(g, 'double')));
                % the matrix with a given Frobenius norm F maximizing the trace is F*eye(d)/sqrt(d),
                % which is of trace F*sqrt(d)
                chiError = sqrt(self.dimension)*self.errorBound; % worst error on the character value
                chiLB = chi - chiError;
                chiUB = chi + chiError;

                % Interval
                %          nTUB      self.dimension
                %  -----------]      X
                %  nontriv. chars    trivial char
                %
                % Case 1 [-]                   -> nontrivial, b = false
                % Case 2 [------]              -> nontrivial, b = false
                % Case 3 [-------------]       -> inconclusive
                % Case 4        [--]           -> inconsistent
                % Case 5        [------]       -> trivial, b = true
                % Case 6                [---]  -> inconsistent
                %
                % ^ plotting [chi-chiError, chi+chiError]
                if chiLB < nonTrivialUB
                    if chiUB < self.dimension
                        % case 1 or 2
                        b = false;
                    else
                        error('Representation is not precise enough to compute the kernel.')
                    end
                elseif chiLB < self.dimension
                    if chiUB < self.dimension
                        % case 4
                        error('Inconsistent character value');
                    else
                        % case 5
                        b = true;
                    end
                else
                    error('Inconsistent character value');
                end
            end
        end

        function K = computeKernel(self)
            assert(isa(self.group, 'replab.FiniteGroup'));
            % TODO error estimation: take in account the uncertainty on computed images
            C = self.group.conjugacyClasses.classes;
            inK = cellfun(@(c) self.inKernel(c.representative), C);
            K = self.group.subgroup(cellfun(@(c) c.representative, C(inK), 'uniform', 0));
            K = self.group.normalClosure(K);
        end

        function f = computeFrobeniusSchurIndicator(self)
            if isa(self.group, 'replab.FiniteGroup')
                f = 0;
                C = self.group.conjugacyClasses.classes;
                n = length(C);
                g2 = cellfun(@(c) self.group.composeN(c.representative, 2), C, 'uniform', 0);
                factor = cellfun(@(c) self.group.order/c.nElements, C, 'uniform', 0);
                if self.canComputeType('cyclotomic')
                    f = replab.cyclotomic.zeros(1, 1);
                    for i = 1:n
                        f = f + trace(self.image(g2{i}, 'cyclotomic'))/replab.cyclotomic.fromVPIs(factor{i});
                    end
                    f = double(f);
                    assert(isreal(f) && round(f) == f);
                    f = round(f);
                else
                    if self.errorBound >= 1
                        error('Error on this representation is too big to compute the Frobenius-Schur indicator');
                    end
                    f = 0;
                    for i = 1:n
                        f = f + trace(self.image(g2{i}, 'double'))/double(factor{i});
                    end
                    f = round(f);
                end
            else
                error('Does not work with continuous groups.');
                % TODO: if this representation has a decomposition, use it
            end

        end

    end

    methods % Representation properties

        function c = conditionNumberEstimate(self)
        % Returns an estimation of the condition number of the change of basis that makes this representation unitary
        %
        % Returns:
        %   double: Condition number of the change of basis matrix
            c = self.cached('conditionNumberEstimate', @() self.computeConditionNumberEstimate);
        end

        function e = errorBound(self)
        % Returns a bound on the approximation error of this representation
        %
        % Returns:
        %   double: There exists an exact representation ``repE`` such that ``||rep(g) - repE(g)||fro <= errorBound``
            e = self.cached('errorBound', @() self.computeErrorBound);
        end

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

        %        function b = isIrreducible(self)
        %        end

        %function d = trivialDimension(self)
        % Returns the dimension of the trivial subrepresentation in this representation
        %
        % Returns:
        %   integer: Dimension
        %
        %end

        function f = frobeniusSchurIndicator(self)
        % Returns the Frobenius-Schur indicator of this representation
        %
        % It corresponds to the value $\iota = \int_{g \in G} tr[\rho_g^2] d \mu$ or
        % $\iota = \frac{1}{|G|} \sum_{g \in G} tr[\rho_g^2]$.
        %
        % Returns:
        %   integer: Value of the indicator
            f = self.cached('frobeniusSchurIndicator', @() self.computeFrobeniusSchurIndicator);
        end

        function K = kernel(self)
        % Returns the kernel of the given representation
        %
        % Only works if `.group` is a finite group.
        %
        % Raises:
        %   An error is this representation is not precise enough to conclude.
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

        %        function b = isDivisionAlgebraCanonical(self)
        %        end


                %   isUnitary (logical): Whether the representation is unitary, mandatory argument
        %   isIrreducible (logical or ``[]``, optional): Whether this representation is irreducible, default ``[]``
        %   trivialDimension (integer or ``[]``, optional): Dimension of the trivial subrepresentation, default ``[]``
        %   frobeniusSchurIndicator % (integer or ``[]``, optional): Exact value of the Frobenius-Schur indicator, default ``[]``
        %   isDivisionAlgebraCanonical % (logical or ``[]``): If the representation is real and irreducible and the Frobenius-Schur indicator is not 1, means that the images encode the complex or quaternion division algebras in the RepLAB canonical form


    end

% $$$     methods (Access = protected)
% $$$
% $$$
% $$$     end
% $$$
% $$$     methods % Computed properties
% $$$
% $$$
% $$$     end

    methods % Derived vector spaces/algebras

% $$$         function e = equivariantTo(self, repR)
% $$$         % Returns the space of equivariant linear maps from this rep to another rep
% $$$         %
% $$$         % The equivariant vector space contains the matrices X such that
% $$$         %
% $$$         % ``X * self.image(g) = repR.image(g) * X``
% $$$         %
% $$$         %
% $$$         % Args:
% $$$         %   repR (`+replab.Rep`): Representation on the target/row space
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Equivariant`: The equivariant vector space
% $$$             e = replab.Equivariant.make(self, repR, '');
% $$$         end
% $$$
% $$$         function c = computeCommutant(self)
% $$$             c = replab.Equivariant.make(self, self, 'commutant');
% $$$         end
% $$$
% $$$         function c = commutant(self)
% $$$         % Returns the commutant of this representation
% $$$         %
% $$$         % This is the algebra of matrices that commute with the representation,
% $$$         % i.e. the vector space isomorphism to the equivariant space from this rep to this rep.
% $$$         %
% $$$         % For any ``g in G``, we have ``rho(g) * X = X * rho(g)``.
% $$$         %
% $$$         % The computation is cached.
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Equivariant`: The commutant algebra represented as an equivariant space
% $$$             c = self.cached('commutant', @() self.computeCommutant);
% $$$         end
% $$$
% $$$         function h = computeHermitianInvariant(self)
% $$$             h = replab.Equivariant.make(self.dual.conj, self, 'hermitian');
% $$$         end
% $$$
% $$$         function h = hermitianInvariant(self)
% $$$         % Returns the Hermitian invariant space of this representation
% $$$         %
% $$$         % This is the space of Hermitian matrices that are invariant under this representation
% $$$         % i.e.
% $$$         %
% $$$         % for any g in G, we have ``X = rho(g) * X * rho(g^-1)'``
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Equivariant`: The space of Hermitian invariant matrices
% $$$             h = self.cached('hermitianInvariant', @() self.computeHermitianInvariant);
% $$$         end
% $$$
% $$$         function t = computeTrivialSpace(self)
% $$$             tRep = self.group.trivialRep(self.field, self.dimension);
% $$$             t = replab.Equivariant.make(tRep, self, 'trivial');
% $$$         end
% $$$
% $$$         function t = trivialSpace(self)
% $$$         % Returns an equivariant space from a trivial representation to this representation
% $$$         %
% $$$         % The trivial representation has the same dimension as this representation
% $$$             t = self.cached('trivialSpace', @() self.computeTrivialSpace);
% $$$         end

    end

% $$$     methods % Irreducible decomposition
% $$$
% $$$         function I = decomposition(self)
% $$$         % Returns the irreducible decomposition of this representation
% $$$         %
% $$$         % Requires this representation to be unitary
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Irreducible`: The irreducible decomposition
% $$$         %
% $$$         % Raises:
% $$$         %   An error is this representation is not unitary.
% $$$             I = self.cached('decomposition', @() self.computeDecomposition);
% $$$         end
% $$$
% $$$         function dec = computeDecomposition(self)
% $$$             dec = replab.irreducible.decomposition(self);
% $$$             if dec.nComponents == 1 && dec.components{1}.multiplicity == 1
% $$$                 assert(~isequal(self.isIrreducible, false));
% $$$                 self.isIrreducible = true;
% $$$                 if isequal(dec.basis, eye(self.dimension))
% $$$                     replab.rep.copyProperties(dec, self);
% $$$                 end
% $$$             end
% $$$         end
% $$$
% $$$     end

    methods % Implementations

        % Str

        function s = headerStr(self)
            p = {};
            if self.isUnitary
                if self.overR
                    p{1,end+1} = 'orthogonal';
                else
                    p{1,end+1} = 'unitary';
                end
            else % ~self.isUnitary
                if self.overR
                    p{1,end+1} = 'nonorthogonal';
                else
                    p{1,end+1} = 'nonunitary';
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
        %   double(\*,\*): The matrix ``M * self.image(inverse(g))``
            gInv = self.group.inverse(g);
            M = full(M * self.image_internal(gInv));
        end

    end

    methods % Morphism composition

% $$$         function res = imap(self, f)
% $$$         % Maps the representation under an isomorphism
% $$$         %
% $$$         % Args:
% $$$         %   f (`.FiniteIsomorphism`): Isomorphism with ``self.group.isSubgroupOf(f.source)``
% $$$         %
% $$$         % Returns:
% $$$         %   `.Rep`: Representation satisfying ``newRep.group.isSubgroup(f.target)``.
% $$$             if self.group.order < f.source.order
% $$$                 f = f.restrictedSource(self.group);
% $$$             end
% $$$             res = replab.rep.CompositionRep(f.inverse, self);
% $$$             res.isUnitary = self.isUnitary;
% $$$             res.trivialDimension = self.trivialDimension;
% $$$             res.isIrreducible = self.isIrreducible;
% $$$             res.frobeniusSchurIndicator = self.frobeniusSchurIndicator;
% $$$             res.isDivisionAlgebraCanonical = self.isDivisionAlgebraCanonical;
% $$$         end

        function res = compose(self, applyFirst)
        % Composition of a representation with a morphism, with the morphism applied first
        %
        % Args:
        %   applyFirst (`.Morphism`): Morphism from a finite group
        %
        % Returns:
        %   `.Rep`: Representation
            res = replab.rep.CompositionRep(applyFirst, self);
        end

    end

    methods % Derived representations

% $$$         function rep1 = contramap(self, morphism)
% $$$         % Returns the representation composed with the given morphism applied first
% $$$         %
% $$$         % Args:
% $$$         %   morphism (`+replab.FiniteMorphism`): Morphism of finite groups such that ``morphism.target == self.group``
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: Representation on the finite group ``morphism.source``
% $$$             assert(self.group == morphism.target);
% $$$             rep1 = replab.Rep.lambda(morphism.source, self.field, self.dimension, @(g) rep.image_internal(morphism.image(g)), @(g) replab.inverseImage_internal(morphism.image(g)));
% $$$         end

% $$$         function rep1 = restrictedTo(self, subgroup)
% $$$         % Returns the restricted representation to the given subgroup
% $$$         %
% $$$         % Args:
% $$$         %   subgroup (`.CompactGroup`): Subgroup to restrict the representation to, must be a subgroup of `.group`
% $$$             rep1 = replab.Rep.lambda(subgroup, self.field, self.dimension, @(g) self.image(g), @(g) self.inverseImage(g));
% $$$         end
% $$$
% $$$         function rep1 = simplify(self)
% $$$         % Returns a representation identical to this, but possibly with its composition simplified
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: A possibly simplified representation
% $$$             rep1 = replab.rep.simplify(self);
% $$$         end
% $$$
% $$$         function complexRep = complexification(self)
% $$$         % Returns the complexification of a real representation
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The complexification of this representation
% $$$         %
% $$$         % Raises:
% $$$         %   An error if this representation is already complex.
% $$$             assert(self.overR, 'Representation should be real to start with');
% $$$             complexRep = replab.rep.ComplexifiedRep(self);
% $$$         end

% $$$         function rep = conj(self)
% $$$         % Returns the complex conjugate representation of this representation
% $$$         %
% $$$         % See https://en.wikipedia.org/wiki/Complex_conjugate_representation
% $$$         %
% $$$         % It obeys ``rep.conj.image(g) = conj(rep.image(g))``
% $$$         %
% $$$         % If this representation is real, it is returned unchanged.
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The complex conjugate of this representation
% $$$             rep = replab.rep.DerivedRep(self, true, false, false);
% $$$         end
% $$$
% $$$         function rep = dual(self)
% $$$         % Returns the dual representation of this representation
% $$$         %
% $$$         % See https://en.wikipedia.org/wiki/Dual_representation
% $$$         %
% $$$         % It obeys ``rep.dual.image(g) = rep.inverseImage(g).'``
% $$$         %
% $$$         % Returns:
% $$$         %   replab.Rep: The dual representation
% $$$             rep = replab.rep.DerivedRep(self, false, true, true);
% $$$         end
% $$$
% $$$         function rep = blkdiag(varargin)
% $$$         % Direct sum of representations
% $$$         %
% $$$         % See `+replab.CompactGroup.directSumRep`
% $$$             self = varargin{1};
% $$$             rep = self.group.directSumRep(self.field, varargin);
% $$$         end
% $$$
% $$$         function rep = kron(varargin)
% $$$         % Tensor product of representations
% $$$         %
% $$$         % See `+replab.CompactGroup.tensorRep`
% $$$             self = varargin{1};
% $$$             rep = self.group.tensorRep(self.field, varargin);
% $$$         end
% $$$
% $$$         function rep = tensorPower(self, n)
% $$$         % Returns a tensor power of this representation
% $$$         %
% $$$         % Args:
% $$$         %   n (integer): Exponent of the tensor power
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The tensor power representation
% $$$             reps = arrayfun(@(i) self, 1:n, 'uniform', 0);
% $$$             rep = self.group.tensorRep(self.field, reps);
% $$$         end
% $$$
% $$$         function rep = directSumOfCopies(self, n)
% $$$         % Returns a direct sum of copies of this representation
% $$$         %
% $$$         % Args:
% $$$         %   n (integer): Number of copies
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The direct sum representation
% $$$             reps = arrayfun(@(i) self, 1:n, 'uniform', 0);
% $$$             rep = self.group.directSum(self.field, reps);
% $$$         end

% $$$         %% Manipulation of representation space
% $$$
% $$$         function res = computeUnitarize(self)
% $$$             if isequal(self.isUnitary, true)
% $$$                 res = replab.SimilarRep.identical(self);
% $$$             else
% $$$                 [A Ainv] = self.unitaryChangeOfBasis;
% $$$                 res = self.similarRep(A, Ainv);
% $$$                 res.isUnitary = true;
% $$$             end
% $$$         end
% $$$
% $$$         function res = unitarize(self)
% $$$         % Returns a unitary representation equivalent to this representation
% $$$         %
% $$$         % The returned representation is of type `.SimilarRep`, from which
% $$$         % the change of basis matrix can be obtained.
% $$$         %
% $$$         % If the representation is already unitary, the returned `.SimilarRep`
% $$$         % has the identity matrix as a change of basis.
% $$$         %
% $$$         % Example:
% $$$         %   >>> S3 = replab.S(3);
% $$$         %   >>> defRep = S3.naturalRep.complexification;
% $$$         %   >>> C = randn(3,3) + 1i * rand(3,3);
% $$$         %   >>> nonUnitaryRep = defRep.subRep(C, inv(C));
% $$$         %   >>> unitaryRep = nonUnitaryRep.unitarize;
% $$$         %   >>> U = unitaryRep.sample;
% $$$         %   >>> norm(U*U' - eye(3)) < 1e-10
% $$$         %      ans =
% $$$         %       logical
% $$$         %       1
% $$$         %
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.SimilarRep`: Unitary similar representation
% $$$             res = self.cached('unitarize', @() self.computeUnitarize);
% $$$         end

% $$$         function [sub1 sub2] = maschke(self, basis1, embedding1)
% $$$         % Given a basis of a subrepresentation, returns two complementary subrepresentations
% $$$         %
% $$$         % Note that we optimize special cases when the representation and/or the basis is
% $$$         % unitary
% $$$         %
% $$$         % Args:
% $$$         %   basis1 (double(dParent,dSub1)): Basis of the first subrepresentation
% $$$         %   embedding1 (double(dChild,dParent), optional): Map from the parent space to the subrepresentation
% $$$         %
% $$$         % Returns
% $$$         % -------
% $$$         %   sub1: `replab.SubRep`
% $$$         %     First subrepresentation
% $$$         %   sub2: `replab.SubRep`
% $$$         %     Second subrepresentation
% $$$             assert(size(basis1, 1) == self.dimension);
% $$$             d1 = size(basis1, 2);
% $$$             if nargin > 2
% $$$                 assert(size(embedding1, 1) == d1);
% $$$                 assert(size(embedding1, 2) == self.dimension);
% $$$             end
% $$$             rest = null(basis1.');
% $$$             if isequal(self.isUnitary, true)
% $$$                 if nargin < 3
% $$$                     BB = basis1'*basis1;
% $$$                     BB = (BB+BB')/2;
% $$$                     embedding1 = BB \ basis1';
% $$$                 end
% $$$                 sub1 = replab.SubRep(self, basis1, embedding1);
% $$$                 % for the second subrepresentation, we take the orthogonal complement from 'null'
% $$$                 sub2 = replab.SubRep(self, rest, rest');
% $$$             else
% $$$                 X = [basis1 rest];
% $$$                 Xinv = inv(X);
% $$$                 P = basis1 * Xinv(1:d1, :);
% $$$                 P1 = self.commutant.project(P);
% $$$                 embedding1 = basis1 \ P1;
% $$$                 P2 = eye(self.dimension) - P1;
% $$$                 basis2 = orth(P2);
% $$$                 embedding2 = basis2 \ P2;
% $$$                 sub1 = replab.SubRep(self, basis1, embedding1);
% $$$                 sub2 = replab.SubRep(self, basis2, embedding2);
% $$$             end
% $$$         end

% $$$         function sub = subRep(self, basis, embedding)
% $$$         % Returns a subrepresentation of this representation
% $$$         %
% $$$         % The subrepresentation is defined by its basis in the parent representation; to compute
% $$$         % images, an embedding map can be provided.
% $$$         %
% $$$         % While the ``basis`` represents in essence a map from the subrepresentation to parent representation,
% $$$         % the embedding map is a map from the parent representation to the subrepresentation.
% $$$         %
% $$$         % The embedding map is not uniquely defined, for example when the subrepresentation contains irreducible
% $$$         % representations that have multiplicites outside the subrepresentation space.
% $$$         %
% $$$         % However, all variants of the embedding map provide identical results when computing images of the
% $$$         % subrepresentation.
% $$$         %
% $$$         % If the embedding is not provided, one is obtained by a trick based on Maschke theorem.
% $$$         % Args:
% $$$         %   basis (double(dParent,dChild), may be sparse): Basis of the subrepresentation
% $$$         %   embedding (double(dChild,dParent), may be sparse, optional): Map from the parent space to the subrepresentation
% $$$         % Returns:
% $$$         %   `+replab.SubRep`: Subrepresentation
% $$$             if nargin < 3
% $$$                 if isequal(self.isUnitary, true)
% $$$                     % optimization for unitary parent representations
% $$$                     embedding = (basis'*basis) \ basis';
% $$$                     % there is an optimization that cannot be made yet:
% $$$                     % if basis'*basis is identity, then embedding = basis'
% $$$                     % works
% $$$                 else
% $$$                     dSub = size(basis, 2);
% $$$                     rest = null(basis.');
% $$$                     X = [basis rest];
% $$$                     Xinv = inv(X);
% $$$                     P = basis * Xinv(1:dSub, :);
% $$$                     P1 = self.commutant.project(P);
% $$$                     embedding = basis \ P1;
% $$$                 end
% $$$             end
% $$$             sub = replab.SubRep(self, basis, embedding);
% $$$         end
% $$$
% $$$         function irreps = splitIntoIrreducibles(self, context)
% $$$         % Decomposes fully the given representation into subrepresentations
% $$$         %
% $$$         % Returns a list of irreducible representations, where trivial subrepresentations
% $$$         % have been identified
% $$$         %
% $$$         % If this representation is irreducible, it will set its `~+replab.Rep.isIrreducible`.
% $$$         %
% $$$         % Args:
% $$$         %   context (`+replab.Context`, optional): Sampling context to use
% $$$         %
% $$$         % Returns:
% $$$         %   cell(1,\*) of `+replab.SubRep`: irreducible subrepresentations
% $$$             if nargin < 2
% $$$                 context = replab.Context.make;
% $$$             end
% $$$             d = self.dimension;
% $$$             start = replab.SubRep.fullSubRep(self);
% $$$             todo = {start};
% $$$             irreps = cell(1, 0);
% $$$             while ~isempty(todo)
% $$$                 h = todo{1};
% $$$                 if isequal(h.isIrreducible, true)
% $$$                     % head of list is irreducible, remove it
% $$$                     irreps{1,end+1} = h;
% $$$                     todo = todo(2:end);
% $$$                 else
% $$$                     res = replab.irreducible.split(h, context);
% $$$                     res1 = cellfun(@(sub) replab.rep.collapse(sub), res, 'uniform', 0);
% $$$                     todo = horzcat(todo(2:end), res1);
% $$$                 end
% $$$             end
% $$$             if nargin < 2
% $$$                 context.close;
% $$$             end
% $$$         end
% $$$
% $$$         function rep1 = similarRep(self, A, Ainv)
% $$$         % Returns a similar representation under a change of basis
% $$$         %
% $$$         % It returns a representation ``rep1`` such that
% $$$         %
% $$$         % ``rep1.image(g) = A * self.image(g) * Ainv``
% $$$         %
% $$$         % Args:
% $$$         %   A (double(\*,\*)): Change of basis matrix
% $$$         %   Ainv (double(\*,\*)): Inverse of the change of basis matrix
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.SimilarRep`: The similar representation
% $$$             rep1 = replab.SimilarRep(self, A, Ainv);
% $$$         end

    end

% $$$     methods (Access = protected)
% $$$
% $$$         function [A Ainv] = unitaryChangeOfBasis(self)
% $$$         % Returns the change of basis to a unitary representation
% $$$         %
% $$$         % Returns ``A`` and ``Ainv`` so that ``A * self.image(g) * Ainv`` is unitary.
% $$$         %
% $$$         % Returns
% $$$         % -------
% $$$         %   A: double(\*,\*)
% $$$         %     Change of basis matrix
% $$$         %   Ainv: double(\*,\*)
% $$$         %     Inverse of change of basis matrix
% $$$             if isequal(self.isUnitary, true)
% $$$                 A = eye(self.dimension);
% $$$                 Ainv = eye(self.dimension);
% $$$             else
% $$$                 X = self.hermitianInvariant.project(eye(self.dimension));
% $$$                 for i = 1:self.dimension
% $$$                     X(i,i) = real(X(i,i));
% $$$                 end
% $$$                 Ainv = chol(X, 'lower');
% $$$                 A = inv(Ainv);
% $$$             end
% $$$         end
% $$$
% $$$     end
% $$$
% $$$     methods (Static)
% $$$
% $$$         function rep = lambda(group, field, dimension, image_internalFun, inverseImage_internalFun)
% $$$         % Creates a non unitary representation from an image function
% $$$         %
% $$$         % Args:
% $$$         %   group (replab.Group): Group represented
% $$$         %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
% $$$         %   dimension (integer): Representation dimension
% $$$         %   image_internalFun (function_handle): Function handle that returns an image matrix given a group element
% $$$         %   inverseImage_internalFun (function_handle): Function handle that returns the inverse of the image
% $$$         %                                               matrix given a group element
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.Rep`: The constructed representation
% $$$             rep = replab.lambda.Rep(group, field, dimension, image_internalFun, inverseImage_internalFun);
% $$$         end
% $$$
% $$$     end

end
