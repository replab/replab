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
%   Either the override the `.image_double_sparse` method if the user implementation provides floating-point
%   images, or both the `.isExact` and `.image_exact` methods.
%
%   The method `.computeErrorBound` must also be overriden (or the value cached at construction).
%
%   This class implements also extra methods about the action of the representation on matrices,
%   methods which can be overloaded for performance.
%
%   If this class has representations as "parts", override the `.decomposeTerm` and `.composeTerm` methods.

    properties (SetAccess = protected)
        group     % (`+replab.CompactGroup`): Group being represented
        field     % ({'R', 'C'}): Vector space defined on real (R) or complex (C) field
        dimension % (integer): Representation dimension
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
        %   isUnitary (logical or ``[]``, optional): Whether the representation is unitary, default ``[]``
        %   isIrreducible (logical or ``[]``, optional): Whether this representation is irreducible, default ``[]``
        %   trivialDimension (integer or ``[]``, optional): Dimension of the trivial subrepresentation, default ``[]``
        %   frobeniusSchurIndicator % (integer or ``[]``, optional): Exact value of the Frobenius-Schur indicator, default ``[]``
        %   isDivisionAlgebraCanonical % (logical or ``[]``): If the representation is real and irreducible and the Frobenius-Schur indicator is not 1, means that the images encode the complex or quaternion division algebras in the RepLAB canonical form
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            args = struct('isUnitary', [], 'isIrreducible', [], 'trivialDimension', [], 'frobeniusSchurIndicator', [], 'isDivisionAlgebraCanonical', []);
            args = replab.util.populateStruct(args, varargin);
            if ~isempty(args.isUnitary)
                self.cache('isUnitary', logical(args.isUnitary), 'error');
            end
            if ~isempty(args.trivialDimension)
                self.cache('trivialDimension', args.trivialDimension, 'error');
            end
            if ~isempty(args.isIrreducible)
                self.cache('isIrreducible', logical(args.isIrreducible), 'error');
            end
            if ~isempty(args.frobeniusSchurIndicator)
                self.cache('frobeniusSchurIndicator', args.frobeniusSchurIndicator, 'error');
            end
            if ~isempty(args.isDivisionAlgebraCanonical)
                self.cache('isDivisionAlgebraCanonical', logical(args.isDivisionAlgebraCanonical), 'error');
            end
        end

    end

    methods (Access = protected) % Protected methods

        function e = computeErrorBound(self)
            error('Abstract');
        end

        function rho = image_double_sparse(self, g)
        % Returns the image of a group element
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %
        % Returns:
        %   double(\*,\*), may be sparse: Element image
            rho = double(self.image_exact(g));
        end

        function rho = image_exact(self, g)
        % Returns the image of a group element
        %
        % Raises:
        %   An error if `.isExact` is false
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %
        % Returns:
        %   cyclotomic(\*,\*): Element image
            error('Exact images not implemented');
        end

    end

    methods % Image computation

        function rho = image(self, g, type)
        % Returns the image of a group element
        %
        % Raises:
        %   An error if ``type`` is ``'exact'`` and `.isExact` is false.
        %
        % Args:
        %   g (element of `.group`): Element being represented
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Image of the given element for this representation
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'double'
                rho = full(self.image_double_sparse(g));
              case 'double/sparse'
                rho = self.image_double_sparse(g);
              case 'exact'
                rho = self.image_exact(g);
              otherwise
                error('Type must be either double, double/sparse or exact');
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
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Image of the inverse of the given element for this representation
            if nargin < 3
                type = 'double';
            end
            if self.cachedOrDefault('isUnitary', false)
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
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
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

    methods (Access = protected) % Simplification

        function c = decomposeTerm(self)
        % Returns the parts of this representation, if any
        %
        % Returns:
        %   cell(1,\*) of `.Rep`: Representations that are part of this representation
            c = cell(1, 0);
        end

        function newRep = composeTerm(self, newParts)
        % Reconstructs the current representation having its parts possibly rewritting
        %
        % Args:
        %   newParts (cell(1,\*) of `.Rep`): New parts, must be equivalent to the current parts up to simplification
        %
        % Returns:
        %   `.Rep`: Reconstructed representation
            newRep = self;
        end

        function copyProperties(self, rep)
        % Copy the cached properties of another representation to this current representation
        %
        % Args:
        %   rep (`.Rep`): Representation equal to the current representation
            if rep.inCache('isUnitary')
                self.cache('isUnitary', rep.isUnitary, '==');
            end
            if rep.inCache('trivialDimension')
                self.cache('trivialDimension', rep.trivialDimension, '==');
            end
            if rep.inCache('isIrreducible')
                self.cache('isIrreducible', rep.isIrreducible, '==');
            end
            if rep.inCache('frobeniusSchurIndicator')
                self.cache('frobeniusSchurIndicator', rep.frobeniusSchurIndicator, '==');
            end
            if rep.inCache('isDivisionAlgebraCanonical')
                self.cache('isDivisionAlgebraCanonical', rep.isDivisionAlgebraCanonical, '==');
            end
            if rep.inCache('kernel')
                self.cache('kernel', rep.kernel, '==');
            end
            if rep.inCache('unitarize')
                u = rep.unitarize;
                u1 = replab.SimilarRep(self, u.A_internal, u.Ainv_internal);
                self.cache('unitarize', u1, 'ignore');
            end
        end

        function [newRep changed] = innermostTerm(self)
        % Performs leftmost innermost tree rewriting
        %
        % Inspired by Algorithm 1 of `<https://doi.org/10.1145/941566.941568>`_ also available at
        % `<https://homepages.cwi.nl/~paulk/publications/TOSEM03.pdf>`_
        %
        % We first make sure the parts of this current term are rewritten, then address the current term.
        %
        % The rewriting steps are implemented as methods with the name ``reduceTerm_XXX``, that return
        % either ``[]`` if the rewrite is unsuccessful, or a rewritten representation if the rewriting step
        % application is successful.
        %
        % Returns
        % -------
        %   newRep:
        %     `.Rep`: Simplified representation
        %   changed:
        %     logical: Whether newRep has been updated
            parts = self.decomposeTerm;
            changed = false;
            for i = 1:length(parts)
                [newPart changed1] = parts{i}.innermostTerm;
                parts{i} = newPart;
                changed = changed | changed1;
            end
            if ~changed
                self1 = self;
            else
                self1 = self.composeTerm(parts);
                self1.copyProperties(self);
            end
            reduced = self1.reduceTerm;
            if isempty(reduced)
                newRep = self1;
            else
                changed = true;
                newRep = reduced;
            end
        end

        function newRep = reduceTerm(self)
        % Attempts a rewriting step
        %
        % Returns:
        %   `.Rep`: Simplified representation or ``[]`` if no rewriting step could be performed
            methodList = replab.compat.methodList(metaclass(self));
            for i = 1:length(methodList)
                M = methodList{i};
                name = M.Name;
                if replab.compat.startsWith(name, 'rewriteTerm_')
                    res = self.(name);
                    if ~isempty(res)
                        fprintf('Applying %s to a representation of dimension %d\n', name, self.dimension);
                        res.copyProperties(self);
                        newRep = res.innermostTerm;
                        return
                    end
                end
            end
            newRep = [];
        end

    end

    methods % Simplification

        function newRep = simplify(self)
        % Returns a representation identical to this, but possibly with its structure simplified
        %
        % Returns:
        %   `+replab.Rep`: A possibly simplified representation
            newRep = self.innermostTerm;
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

        function b = computeIsUnitary(self)
            error('TODO');
        end

        function b = computeIsIrreducible(self)
            error('TODO');
        end

        function d = computeTrivialDimension(self)
            error('TODO');
        end

        function c = computeCommutant(self)
            c = replab.Equivariant.make(self, self, 'commutant');
        end

        function h = computeHermitianInvariant(self)
            h = replab.Equivariant.make(self.dual.conj, self, 'hermitian');
        end

        function t = computeTrivialSpace(self)
            tRep = self.group.trivialRep(self.field, self.dimension);
            t = replab.Equivariant.make(tRep, self, 'trivial');
        end

        function b = computeIsDivisionAlgebraCanonical(self)
            error('TODO');
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

        function d = trivialDimension(self)
        % Returns the dimension of the trivial subrepresentation in this representation
        %
        % Returns:
        %   integer: Dimension
            d = self.cached('trivialDimension', @() self.computeTrivialDimension);
        end

        function f = frobeniusSchurIndicator(self)
        % Returns the Frobenius-Schur indicator of this representation
        %
        % It corresponds to the value $\iota = \int_{g \in G} tr[\rho_g^2] d \mu$ or
        % $\iota = \frac{1}{|G|} \sum_{g \in G} tr[\rho_g^2]$.
        %
        % For real irreducible representations, the Frobenius-Schur indicator can take values:
        % * ``1`` if the representation is of real-type; its complexification is then also irreducible
        % * ``0`` if the representation is of complex-type; it decomposes into two conjugate irreducible
        %   representations over the complex numbers
        % * ``-2`` if the representation is of quaternion-type.
        %
        % Returns:
        %   integer: Value of the indicator
            f = self.cached('frobeniusSchurIndicator', @() self.computeFrobeniusSchurIndicator);
        end

        function b = isIrreducible(self)
        % Returns whether this representation is irreducible over its current field
        %
        % Note that a real representation may become reducible once it is computed over the complex numbers.
        %
        % Returns:
        %   logical: True if this representation is irreducible, false if it has a nontrivial subrepresentation
            b = self.cached('isIrreducible', @() self.computeIsIrreducible);
        end

        function b = isUnitary(self)
        % Returns whether this representation is unitary
        %
        % In the case of approximate representations, we default to ``false``.
        %
        % Returns:
        %   logical: True if this representation is unitary
            b = self.cached('isUnitary', @() self.computeIsUnitary);
        end

        function K = kernel(self)
        % Returns the kernel of the given representation
        %
        % Only works if `.group` is a finite group.
        %
        % Raises:
        %   An error if this representation is not precise enough to compute the kernel.
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

        function b = isDivisionAlgebraCanonical(self)
        % Returns whether the division algebra of this representation is in the RepLAB canonical form
        %
        % This is relevant only for irreducible representations, and is always true for representation over the
        % complex numbers.
        %
        % For irreducible representations over the real numbers, if the Frobenius-Schur indicator is
        %
        % Raises:
        %   An error if this representation is not irreducible.
        %
        % Returns:
        %   logical: Whether the division algebra is canonical
            b = self.cached('isDivisionAlgebraCanonical', @() self.computeIsDivisionAlgebraCanonical);
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

        function t = trivialSpace(self)
        % Returns an equivariant space from a trivial representation to this representation
        %
        % The trivial representation has the same dimension as this representation
            t = self.cached('trivialSpace', @() self.computeTrivialSpace);
        end

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
            if self.cachedOrDefault('isUnitary', false)
                if self.overR
                    p{1,end+1} = 'orthogonal';
                else
                    p{1,end+1} = 'unitary';
                end
            else % ~self.isUnitary
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

        function M = matrixRowAction(self, g, M, type)
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
        %   M (double(\*,\*) or `.cyclotomic`(\*,\*)): Matrix acted upon
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'.
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): The matrix ``self.image(g) * M``
            if nargin < 4
                type = 'double';
            end
            if strcmp(type, 'exact')
                if isa(M, 'double')
                    M = replab.cyclotomic.fromDoubles(M);
                end
                M = self.image(g, 'exact') * M;
            else
                M = self.image(g, 'double/sparse') * M;
                if strcmp(type, 'double')
                    M = full(M);
                end
            end
        end

        function M = matrixColAction(self, g, M, type)
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
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'.
        %
        % Returns:
        %   double(\*,\*): The matrix ``M * self.image(inverse(g))``
            if nargin < 4
                type = 'double';
            end
            gInv = self.group.inverse(g);
            if strcmp(type, 'exact')
                if isa(M, 'double')
                    M = replab.cyclotomic.fromDoubles(M);
                end
                M = M * self.inverseImage(g, 'exact');
            else
                M = M * self.inverseImage(g, 'double/sparse');
                if strcmp(type, 'double')
                    M = full(M);
                end
            end
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
% $$$
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

    end

    methods % Manipulation of representation space

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
        %   >>> nonUnitaryRep = defRep.similarRep(C, inv(C));
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

        function [sub1 sub2] = maschke(self, injection1, projection1)
        % Given an injection into a subrepresentation, returns two complementary subrepresentations
        %
        % Note that we optimize special cases when the representation and/or the injection is unitary
        %
        % Args:
        %   injection (double(D,d) or `cyclotomic`(D,d), may be sparse): Basis / Injection map for the first subrepresentation
        %   projection (double(D,d) or `cyclotomic`(D,d), may be sparse, optional): Projection map for the first subrepresentation
        %
        % Returns
        % -------
        %   sub1: `replab.SubRep`
        %     First subrepresentation
        %   sub2: `replab.SubRep`
        %     Second subrepresentation
            D = self.dimension;
            assert(size(injection1, 1) == D);
            d1 = size(injection1, 2);
            d2 = D - d1;
            if nargin > 2
                assert(size(projection1, 1) == d1);
                assert(size(projection1, 2) == self.dimension);
                sub1 = self.subRep(injection1, 'projection', projection1);
            else
                sub1 = self.subRep(injection1);
            end
            projector1 = sub1.projector;
            projector2 = speye(D) - projector1;
            if self.cachedOrDefault('isUnitary', false) && sub1.cachedOrDefault('isUnitary', false)
                projector2 = (projector2 + projector2')/2;
                [V, D] = replab.numerical.sortedEig(projector2, 'descend', true);
                injection2 = V(:,1:d2);
                projection2 = injection2';
            else
                [U,S,V] = svd(projector2);
                % 1/sqrt(S) U' * X * V * 1/sqrt(S)
                S = diag(S);
                S = S(1:d2);
                F = diag(1./sqrt(S));
                projection2 = F * U(:,1:d2)';
                injection2 = V(:,1:d2) * F;
            end
            sub2 = self.subRep(injection2, 'projection', projection2);
        end

        function sub = subRep(self, injection, varargin)
        % Constructs a subrepresentation of this representation
        %
        % A subrepresentation is specified by an invariant subspace of the current representation.
        %
        % For simplicity, this subspace can be provided as column vectors in a matrix passed as the first
        % argument to this function. This matrix is called ``injection`` for reasons that are clarified below.
        %
        % The user should also provide the keyword argument ``isUnitary``, except when this representation is
        % unitary, and the injection map corresponds to an isometry; then RepLAB deduces that the subrepresentation
        % is unitary as well.
        %
        % More precisely, the subrepresentation is defined by maps between two vector spaces:
        %
        % * the vector space $V$ corresponding to this representation, of dimension $D$,
        % * the vector space $W$ corresponding to the created subrepresentation, of dimension $d$.
        %
        % The two maps are:
        %
        % * the injection $I: W \rightarrow V$ that takes a vector in $W$ and injects it in the parent representation,
        % * the projection $P: V \rightarrow W$ that takes a vector in $V$ and extracts/returns its component in $W$.
        %
        % Correspondingly, the map $I$ is given as a $D x d$ matrix, while the map $P$ is given as a $d x D$ matrix.
        %
        % If only the injection map $I$ is given as an argument, a projection $P$ is computed; it is not necessarily
        % unique. In particular, if the range of $I$ spans irreducible representations with non-trivial multiplicities,
        % the recovered projection is chosen arbitrarily.
        %
        % However, all choices of the projection map provide the same results when computing images of the
        % subrepresentation.
        %
        % If a floating-point approximation $\tilde{I}$ of the injection map is given instead of the exact map $I$,
        % an upper bound on the error can be provided; otherwise, it will be automatically estimated. In that
        % case, the projection map $\tilde{P}$ is also considered as approximate, regardless of whether it is
        % user-provided or computed.
        %
        % The error bound is ``mapErrorBound``, and corresponds to an upper bound on both
        % $|| I \tilde{P} - id ||_F$ and $|| \tilde{I} P - id ||_F$; in the expression, we assume
        % that $I$ and $P$ are the exact injections/projections closest to $\tilde{I}$ and $\tilde{P}$.
        %
        % Args:
        %   injection (double(D,d) or `cyclotomic`(D,d), may be sparse): Basis / Injection map
        %
        % Keyword Args:
        %   projection (double(D,d) or `cyclotomic`(D,d), may be sparse, optional): Projection map
        %   mapErrorBound (double, optional): Upper bound as described above
        %   mapConditionNumberEstimate (double, optional): Upper bound on the condition number of both $P$ and $I$
        %   isUnitary (logical, optional): Whether the resulting representation is unitary, may be omitted
        %
        % Returns:
        %   `+replab.SubRep`: Subrepresentation
            args = struct('projection', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            projection = args.projection;
            isExact = isa(injection, 'replab.cyclotomic') && (isempty(projection) || isa(projection, 'replab.cyclotomic'));
            if isempty(projection)
                if self.cachedOrDefault('isUnitary', false)
                    if isExact
                        % slower because cyclotomic doesn't implement \ or /
                        projection = inv(injection'*injection)*injection';
                    else
                        projection = (injection'*injection)\injection';
                    end
                else
                    % this is a projector on the subspace W
                    % as I*  (inv(I'*I)*I'*I) *inv(I'*I)*I' = I*inv(I'*I)*I'
                    P1 = injection*inv(injection'*injection)*injection';
                    P2 = self.commutant.project(P1);
                    % A\B gives X which is the solution A*X=B
                    % P1 = injection * projection
                    if isExact
                        % slower because cyclotomic doesn't implement \ or /
                        projection = inv(injection'*injection)*injection'*P1;
                    else
                        projection = injection \ P1;
                    end
                end
            end
            sub = replab.SubRep(self, injection, projection, restArgs{:});
        end

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
        function rep1 = similarRep(self, A, varargin)
        % Returns a similar representation under a change of basis
        %
        % It returns a representation ``rep1`` such that
        %
        % ``rep1.image(g) = A * self.image(g) * inv(A)``.
        %
        % Additional keyword arguments will be passed on to the `.Rep` constructor. The argument
        % ``isUnitary`` must be provided, except when the ``inverse`` is provided and both this
        % representation and the change of basis is unitary.
        %
        % Args:
        %   A (double(\*,\*) or `.cyclotomic`(\*,\*)): Change of basis matrix
        %
        % Keywords Args:
        %   inverse (double(\*,\*) or `.cyclotomic`(\*,\*)): Inverse of the change of basis matrix
        %
        % Returns:
        %   `+replab.SimilarRep`: The similar representation
            args = struct('inverse', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            inverse = args.inverse;
            if isempty(inverse)
                inverse = inv(A);
            end
            rep1 = replab.SimilarRep(self, A, inverse, restArgs{:});
        end

    end

    methods (Access = protected)

        function res = computeUnitarize(self)
            if self.cachedOrDefault('isUnitary', false)
                res = replab.SimilarRep.identical(self);
            else
                [A Ainv] = self.unitaryChangeOfBasis;
                res = self.similarRep(A, 'inverse', Ainv, 'isUnitary', true);
            end
        end

        function [A Ainv] = unitaryChangeOfBasis(self)
        % Returns the change of basis to a unitary representation
        %
        % Returns ``A`` and ``Ainv`` so that ``A * self.image(g) * Ainv`` is unitary.
        %
        % Returns
        % -------
        %   A: double(\*,\*), may be sparse
        %     Change of basis matrix
        %   Ainv: double(\*,\*), may be sparse
        %     Inverse of change of basis matrix
            if self.cachedOrDefault('isUnitary', false)
                A = speye(self.dimension);
                Ainv = speye(self.dimension);
            else
                X = self.hermitianInvariant.project(speye(self.dimension), 'double/sparse');
                for i = 1:self.dimension
                    X(i,i) = real(X(i,i));
                end
                Ainv = chol(X, 'lower');
                A = inv(Ainv);
            end
        end

    end
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
