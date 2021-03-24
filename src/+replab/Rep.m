classdef Rep < replab.Obj
% Describes a finite dimensional representation of a compact group
%
% This class has properties that correspond to information that can be computed and cached
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
        isUnitary % (logical): Whether the representation is unitary
        divisionAlgebraName % ('R->C', 'C->R', 'H->C', 'H->R:rep', ''): Name of the division algebra encoded by this representation (see `.DivisionAlgebra`)
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
        %   isUnitary (logical, optional): Whether the representation is unitary, default: false
        %   divisionAlgebraName % ('R->C', 'C->R', 'H->C', 'H->R:rep', '', optional): Name of the division algebra encoding this representation respects, default ``''``
        %   isIrreducible (logical or ``[]``, optional): Whether this representation is irreducible, default ``[]``
        %   trivialDimension (integer or ``[]``, optional): Dimension of the trivial subrepresentation, default ``[]``
        %   frobeniusSchurIndicator (integer or ``[]``, optional): Exact value of the Frobenius-Schur indicator, default ``[]``
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            args = struct('isUnitary', false, 'isIrreducible', [], 'trivialDimension', [], 'frobeniusSchurIndicator', [], 'divisionAlgebraName', '');
            args = replab.util.populateStruct(args, varargin);
            % fill properties
            self.isUnitary = args.isUnitary;
            self.divisionAlgebraName = args.divisionAlgebraName;
            if ~isempty(args.trivialDimension)
                self.cache('trivialDimension', args.trivialDimension, 'error');
            end
            if ~isempty(args.isIrreducible)
                self.cache('isIrreducible', logical(args.isIrreducible), 'error');
            end
            if ~isempty(args.frobeniusSchurIndicator)
                self.cache('frobeniusSchurIndicator', args.frobeniusSchurIndicator, 'error');
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
        % Raises:
        %   An error if ``type`` is ``'exact'`` and `.isExact` is false.
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
            if rep.inCache('trivialDimension')
                self.cache('trivialDimension', rep.trivialDimension, '==');
            end
            if rep.inCache('isIrreducible')
                self.cache('isIrreducible', rep.isIrreducible, '==');
            end
            if rep.inCache('frobeniusSchurIndicator')
                self.cache('frobeniusSchurIndicator', rep.frobeniusSchurIndicator, '==');
            end
            if rep.inCache('divisionAlgebraName')
                self.cache('divisionAlgebraName', rep.divisionAlgebraName, 'isequal');
            end
            if rep.inCache('kernel')
                self.cache('kernel', rep.kernel, '==');
            end
        end

        function [newRep changed] = innermostTerm(self, options)
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
        % Args:
        %   options (struct): Structure of parameters
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
                [newPart changed1] = parts{i}.innermostTerm(options);
                parts{i} = newPart;
                changed = changed | changed1;
            end
            if ~changed
                self1 = self;
            else
                self1 = self.composeTerm(parts);
                self1.copyProperties(self);
            end
            reduced = self1.reduceTerm(options);
            if isempty(reduced)
                newRep = self1;
            else
                changed = true;
                newRep = reduced;
            end
        end

        function newRep = reduceTerm(self, options)
        % Attempts a rewriting step
        %
        % Args:
        %   options (struct): Structure of parameters
        %
        % Returns:
        %   `.Rep`: Simplified representation or ``[]`` if no rewriting step could be performed
            methodList = replab.compat.methodList(metaclass(self));
            for i = 1:length(methodList)
                M = methodList{i};
                name = M.Name;
                if replab.compat.startsWith(name, 'rewriteTerm_')
                    res = self.(name)(options);
                    if ~isempty(res)
                        replab.msg(1, 'Applying %s to a representation of dimension %d\n', name, self.dimension);
                        res.copyProperties(self);
                        newRep = res.innermostTerm(options);
                        return
                    end
                end
            end
            newRep = [];
        end

    end

    methods % Simplification

        function newRep = simplify(self, varargin)
        % Returns a representation identical to this, but possibly with its structure simplified
        %
        % It pushes `.SubRep` to the outer level, and `.DerivedRep` to the inner levels;
        % expands tensor products.
        %
        % Additional keyword arguments can be provided as key-value pairs.
        %
        % Keyword Args:
        %   dense (logical, optional): Whether to allow multiplication of non-sparse matrices
        %   approximate (logical, optional): Whether to allow multiplication of non-exact matrices, implies ``dense = true``
        %
        % Returns:
        %   `+replab.Rep`: A possibly simplified representation
            options = struct('dense', false, 'approximate', false);
            options = replab.util.populateStruct(options, varargin);
            if options.approximate
                options.dense = true;
            end
            newRep = self.innermostTerm(options);
        end

    end

    methods (Access = protected)

        function c = computeConditionNumberEstimate(self)
            if self.isUnitary
                c = 1;
            else
                simRep = self.unitarize;
                c = cond(simRep.injection_internal);
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
        % Computes the Frobenius-Schur indicator
            f = replab.rep.frobeniusSchurIndicator(self);
        end

        function b = computeIsIrreducible(self)
            b = length(self.split) == 1;
        end

        function iso = computeTrivialComponent_exact(self)
        % Computes the maximal trivial subrepresentation of this representation in exact arithmetic
        %
        % Returns:
        %   `.Isotypic`: Subrepresentation as isotypic component
            assert(self.isExact);
            replab.msg(1, '*** Computing trivial component (exact) of representation of dim = %d', self.dimension);
            P2 = self.trivialProjector('exact');
            d = trace(P2);
            assert(iswhole(d));
            d = double(d);
            replab.msg(2, 'Starting exact LU decomposition');
            t = cputime;
            [L, U, p] = lu(P2);
            replab.msg(2, 'Time (LU decomposition): %2.2f s', cputime - t);
            zeroRows = find(arrayfun(@(i) all(U(i,:) == 0), 1:size(U, 1)));
            q = [setdiff(1:size(U, 1), zeroRows) zeroRows];
            U = U(q,:);
            L = L(:,q);
            L = L(p, 1:d);
            U = U(1:d, :);
            % P2 = L * U -> I = L
            I = L;
            P = U;
            if any(any(P*I ~= replab.cyclotomic.eye(d)))
                P = inv(U*L)*P;
            end
            sub = self.subRep(I, 'projection', P);
            iso = replab.Isotypic.fromTrivialSubRep(sub);
        end

        function iso = computeTrivialComponent_double(self)
        % Computes the maximal trivial subrepresentation of this representation in floating-point precision
        %
        % Returns:
        %   `.Isotypic`: Subrepresentation as isotypic component
            replab.msg(1, '*** Computing trivial component (double) of representation of dim = %d', self.dimension);
            [P2, err] = self.trivialProjector('double');
            if err >= 1
                error('Representation is not precise enough to compute the trivial dimension.');
            end
            d = round(trace(P2));
            replab.msg(2, 'Starting rank revealing QR decomposition');
            t = cputime;
            [I, P, p] = replab.numerical.sRRQR_rank(P2, 1.5, d);
            replab.msg(2, 'Time (RRQR decomposition): %2.2f s', cputime - t);
            if self.isUnitary
                sub = self.subRep(I, 'projection', I', 'isUnitary', true);
            else
                P(:,p) = P; % apply permutation
                P = (P*I)\P;
                sub = self.subRep(I, 'projection', P);
            end
            iso = replab.Isotypic.fromTrivialSubRep(sub);
        end

        function proj = computeTrivialProjector_exact(self)
        % Computes the projector into the trivial subrepresentation of this representation
        %
        % Returns:
        %   `.cyclotomic`(\*,\*): Exact projector
            assert(self.isExact);
            replab.msg(1, '*** Computing trivial projector (exact) of representation of dim = %d', self.dimension);
            P = replab.cyclotomic.eye(self.dimension);
            t = cputime;
            P1 = self.trivialRowSpace.project(P, 'exact');
            proj = self.trivialColSpace.project(P1, 'exact');
            replab.msg(2, 'Time (trivial projection): %2.2f s', cputime - t);
        end

        function res = computeTrivialProjector_double(self)
        % Computes the projector into the trivial subrepresentation of this representation
        %
        % Returns:
        %   {double(\*,\*), double}: Approximate projector and estimated error in Frobenius norm
            replab.msg(1, '*** Computing trivial projector (double) of representation of dim = %d', self.dimension);
            if self.inCache('trivialProjector_exact')
                % Shortcut: if an exact trivial projector is available, use it instead
                P = self.trivialProjector('exact');
                [proj, E] = P.doubleApproximation;
                res = {proj, norm(E, 'fro')};
                return
            end
            t = cputime;
            P = speye(self.dimension);
            [P1, E1] = self.trivialRowSpace.project(P, 'double');
            [proj, E2] = self.trivialColSpace.project(P1, 'double');
            err = E1 + E2;
            replab.msg(2, 'Time in projection: %2.2f s', cputime - t);
            res = {proj, err};
        end

        function d = computeTrivialDimension(self)
            c = self.trivialComponent('double');
            d = c.dimension;
        end

    end

    methods % Representation properties

        function p = invariantBlocks(self)
        % Returns a partition of the set ``1:self.dimension`` such that the subsets of coordinates correspond to invariant spaces
        %
        % This method does not guarantee that the finest partition is returned.
        %
        % Example:
        %   >>> G = replab.SignedPermutationGroup.of([1 4 7 2 5 8 3 6 9], [1 -2 -3 -4 5 6 -7 8 9], [1 3 2 4 6 5 -7 -9 -8]);
        %   >>> rep = G.naturalRep;
        %   >>> p = rep.invariantBlocks;
        %   >>> isequal(p.blocks, {[1] [2 3 4 7] [5 6 8 9]})
        %       1
        %
        % Returns:
        %   `.Partition`: Partition of coarse invariant blocks
            if isa(self.group, 'replab.FiniteGroup') && self.group.isTrivial
                blocks = arrayfun(@(i) {i}, 1:self.dimension, 'uniform', 0);
                % each coordinate in its own block
                p = replab.Partition.fromBlocks(blocks);
            elseif isa(self.group, 'replab.FiniteGroup')
                mask = (self.image(self.group.generator(1)) ~= 0);
                for i = 2:self.group.nGenerators
                    mask = mask | (self.image(self.group.generator(i)) ~= 0);
                end
                p = replab.UndirectedGraph.fromAdjacencyMatrix(mask).connectedComponents;
            else
                p = replab.Partition.fromBlocks({1:self.dimension}); % by default returns one block
            end
        end

        function kv = knownProperties(self, keys)
        % Returns the known properties among the set of given keys as a key-value cell array
        %
        % Args:
        %   cell(1,\*) of charstring: Property names
        %
        % Returns:
        %   cell(1,\*): Key-value pairs
            kv = cell(1, 0);
            for i = 1:length(keys)
                if self.inCache(keys{i})
                    kv{1,end+1} = keys{i};
                    kv{1,end+1} = self.cachedOrEmpty(keys{i});
                end
            end
        end

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

        function c = trivialComponent(self, type)
        % Returns the trivial isotypic component present in this representation
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   `.Isotypic`: The trivial isotypic component
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            switch type
              case 'double'
                if self.inCache('trivialComponent_exact')
                    T = self.trivialComponent('exact');
                    sub = self.subRep(T.injection('double'), 'projection', T.projection('double'), 'isUnitary', true);
                    c = replab.Isotypic.fromTrivialSubRep(sub);
                else
                    c = self.cached('trivialComponent_double', @() self.computeTrivialComponent_double);
                end
              case 'exact'
                c = self.cached('trivialComponent_exact', @() self.computeTrivialComponent_exact);
              otherwise
                error('Unknown type');
            end
        end

        function [proj, err] = trivialProjector(self, type)
        % Returns the projector into the trivial component present in this representation
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   proj: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     Projector
        %   err: double
        %     Estimated error in Frobenius norm
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            switch type
              case 'double'
                res = self.cached('trivialProjector_double', @() self.computeTrivialProjector_double);
                proj = res{1};
                err = res{2};
              case 'exact'
                proj = self.cached('trivialProjector_exact', @() self.computeTrivialProjector_exact);
                err = 0;
              otherwise
                error('Unknown type');
            end
        end

        function f = frobeniusSchurIndicator(self)
        % Returns the Frobenius-Schur indicator of this representation
        %
        % It is an integer corresponding to the value $\iota = \int_{g \in G} tr[\rho_g^2] d \mu$ or
        % $\iota = \frac{1}{|G|} \sum_{g \in G} tr[\rho_g^2]$.
        %
        %
        % For real irreducible representations, the Frobenius-Schur indicator can take values:
        % * ``1`` if the representation is of real-type; its complexification is then also irreducible
        % * ``0`` if the representation is of complex-type; it decomposes into two conjugate irreducible
        %   representations over the complex numbers
        % * ``-2`` if the representation is of quaternion-type.
        %
        % For complex irreducible representations, it can take the values:
        % * ``1`` if the representation is of real-type,
        % * ``0`` if the representation is of complex-type, i.e. is not equivalent to its conjugate,
        % * ``-1`` if the representation is of quaternion-type.
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

        function b = isIrreducibleAndCanonical(self)
        % Returns whether this representation is irreducible, and has its division algebra in the canonical encoding
        %
        % This is always true for complex irreps (though we may decide to canonicalize quaternion-type
        % complex representations later). For real irreps, it is true if the Frobenius-Schur indicator is known
        % and matches the encoded `.divisionAlgebraName`.
        %
        % Returns:
        %   logical: True if the representation is irreducible and has its division algebra in the canonical encoding
            if ~self.isIrreducible
                b = false;
            elseif self.overC
                b = true;
            else
                % self.overR
                switch self.frobeniusSchurIndicator
                  case 1
                    b = true; % nothing needs to be done for real-type representations
                  case 0
                    b = strcmp(self.divisionAlgebraName, 'C->R');
                  case -2
                    b = strcmp(self.divisionAlgebraName, 'H->R:rep');
                  otherwise
                    error('Real irreps must have frobeniusSchurIndicator equal to -2,0,1. Here = %d', self.frobeniusSchurIndicator);
                end
            end
        end

        function b = knownIrreducible(self)
        % Returns whether this representation is known to be irreducible; only a true result is significant
        %
        % Returns:
        %   logical: True if `.isIrreducible` is known and is true
            b = self.cachedOrDefault('isIrreducible', false);
        end

        function b = knownReducible(self)
        % Returns whether this representation is known to be reducible; only a true result is significant
        %
        % Returns:
        %   logical: True if `.isIrreducible` is known and is false
            b = ~self.cachedOrDefault('isIrreducible', true);
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

    end

    methods % Equivariant spaces

        function c = antilinearInvariant(self, type)
        % Returns the equivariant space of matrices representing equivariant antilinear maps
        %
        % Let ``F`` be an antilinear map such that ``F(alpha * x) = conj(alpha) * F(x)`` for any scalar ``alpha``
        % and vector ``x`. We describe ``F`` using a matrix ``J`` such that ``F(x) = J * conj(x)``.
        %
        % The following conditions are equivalent.
        %
        % * ``F`` is an equivariant map: ``F(rho.image(g) * x) == rho.image(g) * F(x)``
        % * ``J`` is an equivariant matrix: ``rho.image(g) * J == J * conj(rho.image(g))``
        %
        % The computation is cached.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   `+replab.Equivariant`: The space of equivariant antilinear maps described by an equivariant space of matrices
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['antilinearInvariant_' type], @() self.equivariantFrom(conj(self), 'special', 'antilinear', 'type', type));
        end

        function c = commutant(self, type)
        % Returns the commutant of this representation
        %
        % This is the algebra of matrices that commute with the representation, i.e. the vector space isomorphism to
        % the equivariant space from this rep to this rep.
        %
        % For any ``g in G``, we have ``rho(g) * X = X * rho(g)``.
        %
        % The computation is cached.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   `+replab.Equivariant`: The commutant algebra represented as an equivariant space
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['commutant_' type], @() self.equivariantFrom(self, 'special', 'commutant', 'type', type));
        end

        function e = equivariantFrom(self, repC, varargin)
        % Returns the space of equivariant linear maps from another rep to this rep
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``self.image(g) * X == X * repC.image(g)``
        %
        % Args:
        %   repC (`+replab.Rep`): Representation on the source/column space
        %
        % Keyword Args:
        %   special (charstring or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            e = replab.Equivariant.make(self, repC, varargin{:});
        end

        function e = equivariantTo(self, repR, varargin)
        % Returns the space of equivariant linear maps from this rep to another rep
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repR.image(g) * X == X * self.image(g)``
        %
        % Args:
        %   repR (`+replab.Rep`): Representation on the target/row space
        %
        % Keyword Args:
        %   special (charstring or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            e = replab.Equivariant.make(repR, self, varargin{:});
        end

        function h = sesquilinearInvariant(self, type)
        % Returns the Hermitian invariant space of this representation (deprecated)
        %
        % This is the space of Hermitian matrices that are invariant under this representation
        % i.e.
        %
        % for any g in G, we have ``rho(g) * X = X * rho(g^-1)'``
        %
        % The computation is cached.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant space of Hermitian invariant matrices
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            repR = dual(conj(self));
            h = self.cached(['sesquilinearInvariant_' type], @() self.equivariantTo(repR, 'special', 'sesquilinear', 'type', type));
        end

        function t = trivialRowSpace(self, type)
        % Returns an equivariant space to a trivial representation from this representation
        %
        % The trivial representation has the same dimension as this representation
        %
        % The computation is cached.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant space
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            repR = self.group.trivialRep(self.field, self.dimension);
            t = self.cached(['trivialRowSpace_' type], ...
                            @() self.equivariantTo(repR, 'special', 'trivialRows', 'type', type));
        end

        function t = trivialColSpace(self, type)
        % Returns an equivariant space to this representation from a trivial representation
        %
        % The trivial representation has the same dimension as this representation
        %
        % The computation is cached.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant space
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            repC = self.group.trivialRep(self.field, self.dimension);
            t = self.cached(['trivialColSpace_' type], ...
                            @() self.equivariantFrom(repC, 'special', 'trivialCols', 'type', type));
        end

    end

    methods % Image of maximal torus

        function b = hasTorusImage(self)
        % Returns whether a description of representation of the maximal torus subgroup is available
        %
        % This description is used internally in RepLAB to speed up the group averaging process when computing equivariants.
        %
        % Returns:
        %   logical: True if the call to `.torusImage` succeeds
            b = self.group.hasReconstruction && self.group.maximalTorusDimension == 0; % we have defaults for trivial torus
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
        % Returns a simple description of the representation of the group maximal torus
        %
        % Let:
        % * ``t`` be an element of the maximal torus ``T`` of `.group`,
        % * ``mu`` be the morphism from ``T`` to `.group`,
        % * ``torusMap`` be the integer matrix representing the morphism from ``T`` to ``T1`` (a torus of dimension `.dimension`).
        %
        % We have equality between
        % * ``self.image(mu.imageElement(t))`` and
        % * ``torusInjection * T1.definingRep.image(torusMap * t) * torusProjection``
        %
        % Returns
        % -------
        %   torusMap: integer(d, r)
        %     Exponents, with ``d`` the representation dimension `.dimension` and ``r`` the torus rank
        %   torusInjection: double(d, d), may be sparse
        %     Injection from the torus representation to this representation
        %   torusProjection: double(d, d) may be sparse
        %     Projection from this representation to the torus representation
            assert(self.group.hasReconstruction);
            assert(self.group.maximalTorusDimension == 0);
            torusMap = zeros(self.dimension, 0);
            torusInjection = speye(self.dimension);
            torusProjection = speye(self.dimension);
        end

    end

    methods % Irreducible decomposition

        function I = decomposition(self, type)
        % Returns the irreducible decomposition of this representation
        %
        % The type ``'double/sparse'`` is the same as ``double``.
        %
        % Exact decomposition depends on the availability of an algorithm: for example, we decompose finite groups
        % with an available character table using the Serre projection formulas.
        %
        % Args:
        %   type ('double', 'double/sparse', 'exact', optional): Decomposition type, default: double
        %
        % Returns:
        %   `+replab.Irreducible`: The irreducible decomposition
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            if strcmp(type, 'double')
                I = self.cached('decomposition_double', @() self.computeDecomposition_double);
            else
                I = self.cached('decomposition_exact', @() self.computeDecomposition_exact);
            end
        end

    end

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
                    p{1,end+1} = 'real';
                else
                    p{1,end+1} = 'complex';
                end
            end
            if self.inCache('trivialDimension')
                if self.trivialDimension == self.dimension
                    p{1,end+1} = 'trivial';
                elseif self.trivialDimension == 0
                    p{1,end+1} = 'fully nontrivial';
                else
                    p{1,end+1} = 'nontrivial';
                end
            end
            if self.inCache('isIrreducible')
                if self.isIrreducible
                    p{1,end+1} = 'irreducible';
                else
                    p{1,end+1} = 'reducible';
                end
            end
            if self.inCache('frobeniusSchurIndicator') && self.overR
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
              case 'replab.rep.EncodingRep'
                p{1,end+1} = ['encoded representation (' self.type ')'];
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

    methods (Access = protected) % Implementations

        function M = matrixRowAction_double_sparse(self, g, M)
            M = self.image(g, 'double/sparse') * M;
        end

        function M = matrixColAction_double_sparse(self, g, M)
            M = M * self.inverseImage(g, 'double/sparse');
        end

        function M = matrixRowAction_exact(self, g, M)
            M = self.image(g, 'exact') * M;
        end

        function M = matrixColAction_exact(self, g, M)
            M = M * self.inverseImage(g, 'exact');
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
                    M = replab.cyclotomic(M);
                end
                M = self.matrixRowAction_exact(g, M);
            else
                if isa(M, 'replab.cyclotomic')
                    M = double(M);
                end
                M = self.matrixRowAction_double_sparse(g, M);
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
            if nargin < 4
                type = 'double';
            end
            if strcmp(type, 'exact')
                if isa(M, 'double')
                    M = replab.cyclotomic(M);
                end
                M = self.matrixColAction_exact(g, M);
            else
                if isa(M, 'replab.cyclotomic')
                    M = double(M);
                end
                M = self.matrixColAction_double_sparse(g, M);
                if strcmp(type, 'double')
                    M = full(M);
                end
            end
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
        end

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

        function rep1 = contramap(self, morphism)
        % Returns the representation composed with the given morphism applied first
        %
        % Args:
        %   morphism (`+replab.FiniteMorphism`): Morphism of finite groups such that ``morphism.target == self.group``
        %
        % Returns:
        %   `+replab.Rep`: Representation on the finite group ``morphism.source``
            assert(self.group == morphism.target);
            rep1 = replab.Rep.CompositionRep(morphism, self);
        end

        function rep1 = restrictedTo(self, subgroup)
        % Returns the restricted representation to the given subgroup
        %
        % Args:
        %   subgroup (`.CompactGroup`): Subgroup to restrict the representation to, must be a subgroup of `.group`
            rep1 = replab.Rep.lambda(subgroup, self.field, self.dimension, @(g) self.image(g), @(g) self.inverseImage(g));
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
            complexRep = replab.rep.EncodedRep(self, 'R^d -> C^d');
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

        function sub = symmetricSquare(self)
        % Returns the symmetric square of this representation
        %
        % Let $V$ be the vector space corresponding to this representation, and $V \otimes V$ be
        % its second tensor power. Let $T$ be the linear map such that $T(v_1 \otimes v_2) = v_2 \otimes v_1$.
        %
        % Then the symmetric square $S$ is the subspace $S subset V \otimes V$ such that $T(s) = s$ for
        % all $s \in S$.
        %
        % Returns:
        %   `.SubRep`: A subrepresentation of the second tensor power of this representation
            d = self.dimension;
            D = d*(d+1)/2;
            ind = 1;
            injection = zeros(d*d, D);
            for i = 1:d
                injection(i+(i-1)*d, ind) = 1;
                ind = ind + 1;
                for j = i+1:d
                    injection(j+(i-1)*d, ind) = 1;
                    injection(i+(j-1)*d, ind) = 1;
                    ind = ind + 1;
                end
            end
            sub = self.tensorPower(2).subRep(replab.cyclotomic(injection));
        end

        function sub = alternatingSquare(self)
        % Returns the alternating square of this representation
        %
        % Let $V$ be the vector space corresponding to this representation, and $V \otimes V$ be
        % its second tensor power. Let $T$ be the linear map such that $T(v_1 \otimes v_2) = v_2 \otimes v_1$.
        %
        % Then the alternating square $A$ is the subspace $A subset V \otimes V$ such that $T(a) = -a$ for
        % all $a \in A$.
        %
        % Returns:
        %   `.SubRep`: A subrepresentation of the second tensor power of this representation
            d = self.dimension;
            D = d*(d-1)/2;
            ind = 1;
            injection = zeros(d*d, D);
            for i = 1:d
                for j = i+1:d
                    injection(j+(i-1)*d, ind) = 1;
                    injection(i+(j-1)*d, ind) = -1;
                    ind = ind + 1;
                end
            end
            sub = self.tensorPower(2).subRep(replab.cyclotomic(injection));
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

        function res = identifyIrrep(self)
        % Identifies the irrep type and returns an equivalent representation with explicit structure
        %
        % This representation must be irreducible.
        %
        % Example:
        %   >>> Q = replab.QuaternionGroup();
        %   >>> rep = Q.naturalRep.complexification;
        %   >>> irreps = rep.split;
        %   >>> assert(length(irreps) == 2); % splits in two equivalent irreps
        %   >>> sub = irreps{1};
        %   >>> sub1 = sub.identifyIrrep;
        %   >>> strcmp(sub1.divisionAlgebraName, 'H->C')
        %       1
        %
        % Returns:
        %   `+replab.SubRep`: Similar representation with the canonical division algebra encoding
            assert(self.isIrreducible, 'This method only applies to irreducible representations');
            if self.overR
                irreps = replab.irreducible.identifyIrrepsInParent(replab.SubRep.identical(self), self.commutant.sample);
                assert(length(irreps) == 1);
                res = irreps{1};
            else
                res = replab.irreducible.identifyComplexIrrepInParent(replab.SubRep.identical(self), self.antilinearInvariant.sample);
            end
        end

        function res = unitarize(self)
        % Returns a unitary representation equivalent to this representation
        %
        % The returned representation is of type `.SubRep`, from which the change of basis matrix
        % can be obtained.
        %
        % If the representation is already unitary, the returned `.SubRep` has the identity matrix as
        % injection/projection maps;
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> defRep = S3.naturalRep.complexification;
        %   >>> C = randn(3,3) + 1i * rand(3,3);
        %   >>> nonUnitaryRep = defRep.similarRep(C);
        %   >>> unitaryRep = nonUnitaryRep.unitarize;
        %   >>> U = unitaryRep.sample;
        %   >>> norm(U*U' - eye(3)) < 1e-10
        %      ans =
        %       logical
        %       1
        %
        % Returns:
        %   `+replab.SubRep`: Unitary similar representation
            res = self.cached('unitarize', @() self.computeUnitarize);
        end

        function sub2 = maschke(self, sub1)
        % Given a subrepresentation, returns the complementary subrepresentation
        %
        % Note that we optimize special cases when the representation and/or the injection is unitary
        %
        % Args:
        %   sub1 (`.SubRep`): Subrepresentation of this representation
        %
        % Returns:
        %   `.SubRep`: Complementary subrepresentation
            D = self.dimension;
            d1 = sub1.dimension;
            d2 = D - d1;
            proj1 = sub1.projector;
            proj2 = speye(D) - proj1;
            [I, P, p] = replab.numerical.sRRQR_rank(proj2, 1.5, d2);
            if self.isUnitary && sub1.isUnitary
                sub2 = self.subRep(I, 'projection', I', 'isUnitary', true);
            else
                P(:,p) = P;
                P = (P*I)\P;
                sub2 = self.subRep(I, 'projection', P);
            end
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
        %   largeScale (logical, optional): Whether to use the large-scale version of the algorithm, default automatic selection
        %   tolerances (`.Tolerances`): Termination criteria
        %   nSamples (integer, optional): Number of samples to use in the large-scale version of the algorithm, default ``5``
        %
        % Returns:
        %   `+replab.SubRep`: Subrepresentation
            args = struct('projection', [], 'largeScale', self.dimension > 1000, 'tolerances', replab.rep.Tolerances, 'nSamples', 5);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            projection = args.projection;
            % convert cell-encoded cyclotomic matrices, when relevant
            if iscell(injection)
                injection = replab.cyclotomic(injection);
            end
            if iscell(projection)
                projection = replab.cyclotomic(projection);
            end
            isExact = isa(injection, 'replab.cyclotomic') && (isempty(projection) || isa(projection, 'replab.cyclotomic'));
            if isempty(projection)
                if isExact
                    if self.isUnitary
                        % slower because cyclotomic doesn't implement \ or /
                        projection = inv(injection'*injection)*injection';
                    else
                        % this is a projector on the subspace W
                        % as I*  (inv(I'*I)*I'*I) *inv(I'*I)*I' = I*inv(I'*I)*I'
                        P1 = injection*inv(injection'*injection)*injection';
                        P2 = self.commutant.project(P1);
                        % A\B gives X which is the solution A*X=B
                        % P1 = injection * projection
                        % slower because cyclotomic doesn't implement \ or /
                        projection = inv(injection'*injection)*injection'*P2;
                    end
                else % non exact
                    if self.isUnitary
                        projection = (injection'*injection)\injection';
                    else
                        if args.largeScale
                            projection = replab.rep.findProjection_largeScale(self, injection, args.nSamples, args.tolerances, [], []);
                        else
                            P1 = injection*inv(injection'*injection)*injection';
                            P2 = self.commutant.project(P1);
                            projection = injection \ P2;
                        end
                    end
                end
            end
            sub = replab.SubRep(self, injection, projection, restArgs{:});
        end

        function sub = blockSubRep(self, block)
        % Constructs a subrepresentation from a subset of Euclidean coordinates
        %
        % The ``block`` argument can be obtained from `.invariantBlocks`, and this method
        % constructs the sparse injection and projection maps with which to call `.subRep`.
        %
        % Args:
        %   block (integer(1,\*)): Euclidean coordinates of the block
        %
        % Returns:
        %   `.SubRep`: Subrepresentation
            d = length(block);
            injection = sparse(block, 1:d, ones(1, d), self.dimension, d);
            projection = injection';
            sub = self.subRep(injection, 'projection', projection, 'isUnitary', self.isUnitary);
        end

        function irreps = split(self, varargin)
        % Splits this representation into irreducible subrepresentations
        %
        % Keyword Args:
        %   forceNonUnitaryAlgorithms (logical, optional): Whether to force the use of algorithms for not necessarily unitary representations, default: false
            args = struct('forceNonUnitaryAlgorithms', false);
            args = replab.util.populateStruct(args, varargin);
            sample1 = self.commutant.sample;
            sample2 = self.commutant.sample;
            [trivial, nonTrivialIrreps] = replab.irreducible.split(self, sample1, sample2, args.forceNonUnitaryAlgorithms);
            tiso = replab.Isotypic.fromTrivialSubRep(trivial);
            irreps = horzcat(tiso.irreps, nonTrivialIrreps);
        end

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
        %   `+replab.SubRep`: The similar representation
            args = struct('inverse', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            inverse = args.inverse;
            if isempty(inverse)
                inverse = inv(A);
            end
            rep1 = self.subRep(inverse, 'projection', A, restArgs{:});
        end

    end

    methods (Access = protected)

        function dec = computeDecomposition_exact(self)
        % Computes the representation decomposition, exact variant, only available for finite groups
            assert(isa(self.group, 'replab.FiniteGroup'));
            ct = self.group.characterTable;
            assert(ct.hasIrreps, 'All irreps must have explicit constructions available');
            irreps = ct.irreps;
            assert(all(cellfun(@(ir) ir.isExact, irreps)), 'All irreps must be available in exact form');
            if self.overC
                % easy case, use Serre's projection formulas
                components = cellfun(@(irrep) replab.irreducible.isotypicComponent(self, irrep, 'exact'), irreps, 'uniform', 0);
                dec = replab.Irreducible(self, components);
            else
                % over the reals, we first compute the complex decomposition
                assert(all(cellfun(@(ir) ir.frobeniusSchurIndicator == 1, irreps)));
                for i = 1:self.group.nGenerators
                    g = self.group.generator(i);
                    for j = 1:length(irreps)
                        img = irreps{j}.image(g, 'exact');
                        assert(isreal(img), 'All representations must have real images');
                    end
                end
                irreps = cellfun(@(ir) replab.rep.EncodedRep(ir, 'C^d -> R^d'), irreps, 'uniform', 0);
            end
        end

        function dec = computeDecomposition_double(self)
        % Computes the representation decomposition
        %
        % First it splits the representation into irreducibles, before recognizing which
        % irreducible representations are part of the same isotypic component.
            sample1 = self.commutant.sample;
            sample2 = self.commutant.sample;
            forceNonUnitaryAlgorithms = false;
            [trivial, nonTrivialIrreps] = replab.irreducible.split(self, sample1, sample2, forceNonUnitaryAlgorithms);
            tiso = replab.Isotypic.fromTrivialSubRep(trivial);
            [iso, zeroErrors, nonZeroErrors] = replab.irreducible.findIsotypic(self, nonTrivialIrreps, sample2);
            if tiso.dimension > 0
                dec = replab.Irreducible(self, horzcat({tiso}, iso));
            else
                dec = replab.Irreducible(self, iso);
            end
        end

        function res = computeUnitarize(self)
            if self.isUnitary
                res = replab.SubRep.identical(self);
            else
                [A, Ainv] = self.unitaryChangeOfBasis;
                res = self.subRep(Ainv, 'projection', A, 'isUnitary', true);
            end
        end

        function [A, Ainv] = unitaryChangeOfBasis(self)
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
            if self.isUnitary
                A = speye(self.dimension);
                Ainv = speye(self.dimension);
            else
                X = self.sesquilinearInvariant.project(eye(self.dimension), 'double');
                X = (X + X')/2;
                A = chol(X);
                Ainv = inv(A);
            end
        end

    end

end
