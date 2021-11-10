classdef SubRep < replab.Rep
% Describes a subrepresentation of a finite dimensional representation
%
% A subrepresentation corresponds to a subspace $W$ of a parent representation associated
% with the vector space $V$.
%
% For simplicity, we say that the subspace $W$ is specified using a set of basis vectors, amalgamated
% as column vectors in a matrix `.injection`.
%
% More precisely, the subrepresentation is defined by two maps, an injection and a projection.
%
% * The injection map $I: W \rightarrow V$ expresses a vector in $V$ from a vector in $W$.
%
% * The projection map $P: V \rightarrow W$ extracts the $W$ part from a vector in the parent space $V$.
%
% Those maps obey $P I = id$.
%
% If one wants an operator $\Pi: V \rightarrow V$ acting as a projector ($\Pi^2 = \Pi$), then
% one simply defines $\Pi = I P$.
%
% Note that often, only approximations $\tilde{P}$ and $\tilde{I}$ are known by RepLAB. We thus
% store a measure of the error associated with those maps in `.projectorErrorBound`.
%
% Note:
%   The terminology of "injection" and "projection" linear maps was inspired by https://arxiv.org/pdf/0812.4969.pdf

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation of dimension $D$
        injection_internal % (double(D,d) or `.cyclotomic`(D,d), may be sparse): Injection map
        projection_internal % (double(d,D) or `.cyclotomic`(d,D), may be sparse): Projection map
        mapsAreAdjoint % (logical): True if `.parent` is unitary (so the Hermitian adjoint makes sense)  and `.injection` is the conjugate transpose of `.projection`
        isSimilarRep % (logical): True if this `.SubRep` encodes a similarity transformation
    end

    methods

        function self = SubRep(parent, injection_internal, projection_internal, varargin)
        % Constructs a subrepresentation of a parent representation
        %
        % Additional keyword arguments are passed to the `.Rep` constructor.
        %
        % If both `.parent` is unitary, and `.injection_internal` is adjoint to `.projection_internal`, then
        % the `.isUnitary` argument may be omitted as it can be deduced to be true.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of dimension $D$
        %   injection_internal (double(D,d) or `.cyclotomic`(D,d), may be sparse): Injection map $I$
        %   projection_internal (double(d,D) or `.cyclotomic`(d,D), may be sparse): Projection map $P$
        %
        % Keyword Args:
        %   isUnitary (logical, optional): Whether the resulting representation is unitary, may be omitted (see above)
        %   projectorErrorBound (double, optional): Upper bound on || I P - \tilde{I} \tilde{P} ||_F
        %   injectionConditionNumberEstimate (double, optional): Upper bound of the condition number of $\tilde{I}$ (and thus $\tilde{P}$)
            D = size(injection_internal, 1);
            d = size(injection_internal, 2);
            isSimilarRep = (d == D);
            assert(size(projection_internal, 1) == d);
            assert(size(projection_internal, 2) == D);
            assert(parent.dimension == D, 'Incorrect dimension');
            if parent.overR
                assert(isreal(injection_internal) && isreal(projection_internal), 'A real Rep can only have real subrepresentations.');
            end
            args = struct('injectionConditionNumberEstimate', [], 'projectorErrorBound', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            mapsAreAdjoint = parent.isUnitary && all(all(injection_internal == projection_internal'));
            if mapsAreAdjoint
                restArgs = replab.util.keyValuePairsUpdate(restArgs, 'isUnitary', true);
            end
            if parent.inCache('trivialDimension')
                if isSimilarRep
                    [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', parent.trivialDimension);
                    assert(~exists || oldValue == parent.trivialDimension);
                elseif parent.trivialDimension == 0
                    [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', 0);
                    assert(~exists || oldValue == 0);
                end
                if parent.trivialDimension == parent.dimension
                    [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', d);
                    assert(~exists || oldValue == d);
                end
            end
            if isSimilarRep
                if parent.inCache('isIrreducible')
                    [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'isIrreducible', parent.isIrreducible);
                    assert(~exists || oldValue == parent.isIrreducible);
                end
                if parent.inCache('frobeniusSchurIndicator')
                    [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'frobeniusSchurIndicator', parent.frobeniusSchurIndicator);
                    assert(~exists || oldValue == parent.frobeniusSchurIndicator);
                end
            end
            self@replab.Rep(parent.group, parent.field, d, restArgs{:});
            self.parent = parent;
            self.injection_internal = injection_internal;
            self.projection_internal = projection_internal;
            self.mapsAreAdjoint = mapsAreAdjoint;
            self.isSimilarRep = isSimilarRep;
            if ~isempty(args.projectorErrorBound)
                self.cache('projectorErrorBound', args.projectorErrorBound, 'error');
            end
            if ~isempty(args.injectionConditionNumberEstimate)
                self.cache('injectionConditionNumberEstimate', args.injectionConditionNumberEstimate, 'error');
            end
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_isIdentity(self, options)
            if self.mapsAreIntegerValued && all(all(self.injection_internal == speye(self.dimension))) && all(all(self.projection_internal == speye(self.dimension)))
                res = self.parent;
            else
                res = [];
            end
        end

        function res = rewriteTerm_SubRepOfSubRep(self, options)
            if isa(self.parent, 'replab.SubRep')
                if self.mapsAreIntegerValued || self.parent.mapsAreIntegerValued || ...
                        (options.dense && (options.approximate || (self.mapsAreExact && self.parent.mapsAreExact)))
                    newI = self.parent.injection_internal * self.injection_internal;
                    newP = self.projection_internal & self.parent.projection_internal;
                    res = replab.SubRep(self.parent.parent, newI, newP);
                    return
                end
            end
            res = [];
        end

    end

    methods

        function sub1 = withUpdatedMaps(self, injection, projection, varargin)
        % Returns a copy of this subrepresentation with updated injection and projection maps
        %
        % The following properties are copied from the original subrepresentation:
        %
        % * `.isIrreducible`
        % * `.frobeniusSchurIndicator`
        % * `.kernel`
        % * `.trivialDimension`
        %
        % Additional keywords arguments are passed to the `.subRep` method of `.parent`.
        %
        % Args:
        %   injection (double(D,d) or `.cyclotomic`(D,d), may be sparse): Injection map
        %   projection (double(d,D) or `.cyclotomic`(d,D), may be sparse): Projection map
        %
        % Returns:
        %   `.SubRep`: The updated subrepresentation
            sub1 = self.parent.subRep(injection, 'projection', projection, varargin{:});
            if self.inCache('isIrreducible')
                sub1.cache('isIrreducible', self.isIrreducible, '==');
            end
            if self.inCache('frobeniusSchurIndicator')
                sub1.cache('frobeniusSchurIndicator', self.frobeniusSchurIndicator, '==');
            end
            if self.inCache('kernel')
                sub1.cache('kernel', self.kernel, '==');
            end
            if self.inCache('trivialDimension')
                sub1.cache('trivialDimension', self.trivialDimension, '==');
            end
        end

        function sub1 = withNoise(self, injectionMapNoise, projectionMapNoise)
        % Adds Gaussian noise to the injection/projection maps
        %
        % This method has two call conventions.
        %
        % If ``sub.withNoise(sigma)`` is called, then the injection map is changed to
        % ``injection1 = injection + Delta``, where ``Delta`` is a matrix with normally distributed
        % entries of standard deviation ``sigma``. ``Delta`` is real (resp. complex) if the representation
        % is real (resp. complex). Then the original `.projection` map is ignored and ``projection1`` is recovered
        % as in `+replab.Rep.subRep`.
        %
        % If ``sub.withNoise(sigmaI, sigmaP)`` is called, then noise of magnitude ``sigmaI`` is added to the injection map
        % and noise of magnitude ``sigmaP`` is added to the projection map. Then the new noisy maps are corrected
        % so that ``injection1 * projection1`` is still a projector.
        %
        % Args:
        %   injectionMapNoise (double): Magnitude of the noise to add to the injection map
        %   projectionMapNoise (double, optional): Magnitude of the noise to add to projection map
        %
        % Returns:
        %   `.SubRep`: A noisy subrepresentation
            d = self.dimension;
            D = self.parent.dimension;
            switch nargin
              case 2
                if self.overR
                    injection1 = self.injection('double') + randn(D, d)*injectionMapNoise;
                else
                    injection1 = self.injection('double') + (randn(D, d) + 1i*rand(D, d))*injectionMapNoise/sqrt(2);
                end
                sub1 = self.parent.subRep(injection1);
              case 3
                if self.overR
                    injection1 = self.injection('double') + randn(D, d)*injectionMapNoise;
                    projection1 = self.projection('double') + randn(d, D)*projectionMapNoise;
                else
                    injection1 = self.injection('double') + (randn(D, d) + 1i*rand(D, d))*injectionMapNoise/sqrt(2);
                    projection1 = self.projection('double') + (randn(d, D) + 1i*rand(d, D))*projectionMapNoise/sqrt(2);
                end
                projection1 = (projection1 * injection1) \ projection1;
                sub1 = self.parent.subRep(injection1, 'projection', projection1, 'divisionAlgebraName', self.divisionAlgebraName);
              otherwise
                error('Wrong calling convention');
            end
        end

        function b = injectionIsExact(self)
        % Returns whether the injection map is written either with cyclotomics or integer coefficients
        %
        % Returns:
        %   logical: True if that is the case
            I = self.injection_internal;
            b = (isa(I, 'cyclotomic') || replab.numerical.isExact(I));
        end

        function b = mapsAreExact(self)
        % Returns whether both the injection map and the projection map are written either with cyclotomics or integer coefficients
        %
        % Returns:
        %   logical: True if both `.injection` and `.projection` are cyclotomic matrices or have Gaussian integer entries
            I = self.injection_internal;
            P = self.projection_internal;
            b = (isa(I, 'cyclotomic') || replab.numerical.isExact(I)) && ...
                (isa(P, 'cyclotomic') || replab.numerical.isExact(P));
        end

        function b = mapsAreIntegerValued(self)
        % Returns whether both the injection map and the projection map are expressed with Gaussian integer coefficients
        %
        % Returns:
        %   logical: True if both `.injection` and `.projection` have Gaussian integer entries
            I = self.injection_internal;
            P = self.projection_internal;
            b = self.mapsAreExact && all(all(I == round(double(I)))) && all(all(P == round(double(P))));
        end

        function e = projectorErrorBound(self)
        % Returns an upper bound on the quality of the approximation of the injection/projection maps
        %
        % Let $I$ and $P$ be the respective exact injection and projection maps. Let $\tilde{I}$ and
        % $\tilde{P}$ be their approximate counterparts.
        %
        % This returns an upper bound on the Frobenius norm of $I P - \tilde{I} \tilde{P}$.
        %
        % Returns:
        %   double: Upper bound on the Frobenius distance of the approximate projector to the commutant space
            e = self.cached('projectorErrorBound', @() self.computeProjectorErrorBound);
        end

        function e = biorthogonalityErrorBound(self)
        % Returns an upper bound on the biorthogonality of the injection/projection maps
        %
        % This returns an upper bound on ``norm(self.projection*self.injection - eye(d), 'fro')``, where
        % ``d == self.dimension``.
            e = self.cached('biorthogonalityErrorBound', @() norm(self.projection*self.injection - eye(self.dimension), 'fro'));
        end

        function c = injectionConditionNumberEstimate(self)
        % Returns an upper bound on the condition number of both the injection and the projection map
        %
        % Returns:
        %   double: Upper bound on the condition number
            c = self.cached('injectionConditionNumberEstimate', @() self.computeInjectionConditionNumberEstimate);
        end

        function mat = injection(self, type)
        % Returns the injection from the subrepresentation to the parent representation
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): The injection map
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            mat = replab.numerical.convert(self.injection_internal, type);
        end

        function mat = basis(self, indices, type)
        % Returns (part of) the basis of the subrepresentation in the parent representation
        %
        % It is an alias for ``I(:, indices)`` where ``I = self.injection(type)``.
        %
        % Args:
        %   indices (integer(1,\*) or ``[]``): Indices of the basis elements to return
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): The basis of the subrepresentation given as column vectors
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if nargin < 2 || isempty(indices)
                indices = 1:self.dimension;
            end
            mat = self.injection(type);
            mat = mat(:, indices);
        end

        function mat = projection(self, type)
        % Returns the projection from the parent representation to the subrepresentation
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): The projection map
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            mat = replab.numerical.convert(self.projection_internal, type);
        end

        function mat = projector(self, type)
        % Returns the projector on this subrepresentation
        %
        % Note that the projector is in general not uniquely defined by the injection map;
        % it is defined by the injection map only when the subrepresentation is a direct sum
        % of isotypic components.
        %
        % Returns:
        %   double(\*,\*): Projector matrix on the subrepresentation
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            mat = self.injection(type) * self.projection(type);
        end

        function [res, better] = nice(self)
        % Returns a representation similar to the current subrepresentation, with a nicer basis (EXPERIMENTAL)
        %
        % The "niceness" of the basis is implementation dependent. As of the first implementation
        % of this feature, RepLAB tries to make the basis real, and then with small integer
        % coefficients.
        %
        % The returned subrepresentation is not necessarily unitary.
        %
        % In the case no improvement could be made, the original subrepresentation is returned.
        %
        % Returns
        % -------
        %   res: `.SubRep`
        %      A subrepresentation of ``self.parent``
        %   better: logical
        %      Whether a nicer subrepresentation has been found
            res = replab.rep.niceSubRep(self);
            better = ~isempty(res);
            if ~better
                res = self;
            end
        end

        function sub1 = refine(self, varargin)
        % Refines this subrepresentation
        %
        % Assumes that the subrepresentation is already close to an exact subrepresentation, and refines its subspace
        % by an iterative procedure applied on its `.injection` and `.projection` maps.
        %
        % Keyword Args:
        %   tolerances (`+replab.+rep.Tolerances`): Termination criteria
        %   largeScale (logical, optional): Whether to use the large-scale version of the algorithm, default automatic choice
        %   nSamples (integer, optional): Number of samples to use in the large-scale version of the algorithm, default ``5``
        %   nInnerIterations (integer, optional): Number of inner iterations in the medium-scale version of the algorithm, default ``3``
        %   injectionBiortho (double(\*,\*), may be sparse): Injection map of known multiplicity space to remove from this subrepresentation
        %   projectionBiortho (double(\*,\*), may be sparse): Projection map of known multiplicity space to remove from this subrepresentation
        %
        % Returns
        % -------
        %   sub1: `.SubRep`
        %     Subrepresentation with refined subspace (injection/projection maps)
            gen = replab.rep.GenSubRep.fromSubRep(self);
            gen1 = gen.refine(varargin{:});
            res = gen1.toSubRep;
            sub1 = self.withUpdatedMaps(res.injection, res.projection, 'isUnitary', res.isUnitary, 'divisionAlgebraName', self.divisionAlgebraName);
        end

        function res = collapse(self)
        % Simplifies a SubRep of a SubRep
        %
        % Raises:
        %   An error if `.parent` is not of type `.SubRep`
        %
        % Returns:
        %   `.SubRep`: A subrepresentation of ``.parent.parent``
            parent = self.parent;
            switch class(parent)
              case 'replab.SubRep'
                newI_internal = parent.injection_internal * self.injection_internal;
                newP_internal = self.projection_internal * parent.projection_internal;
                res = parent.parent.subRep(newI_internal, 'projection', newP_internal);
                res.copyProperties(self);
              otherwise
                error('Not supported');
            end
        end

    end

    methods % Equivariant spaces

        function E = subEquivariantFrom(self, repC, varargin)
        % Returns the space of equivariant linear maps from another subrepresentation to this subrepresentation
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * repC.image(g) = self.image(g) * X``
        %
        % Example:
%         %   >>> S3 = replab.S(3);
%         %   >>> rep = S3.naturalRep;
%         %   >>> Xrep = randn(3, 3);
%         %   >>> triv = rep.subRep([1;1;1]);
%         %   >>> std = rep.maschke(triv);
%         %   >>> E = triv.subEquivariantFrom(std);
%         %   >>> Xsub = E.projectFromParent(Xrep);
%         %   >>> isequal(size(Xsub), [1 2]) % map to trivial (dim=1) from std (dim=2)
%         %       1
%         %   >>> norm(Xsub) <= 1e-15
%         %       1
        %
        % Args:
        %   repC (`+replab.SubRep`): Subrepresentation on the source/column space
        %
        % Keyword Args:
        %   special (charstring, optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %   parent (`.Equivariant`, optional): Equivariant space of ``self.parent`` and ``repC.parent``, default: ``[]``
        %
        % Returns:
        %   `+replab.SubEquivariant`: The equivariant vector space
            E = replab.SubEquivariant.make(self, repC, varargin{:});
        end

        function E = subEquivariantTo(self, repR, varargin)
        % Returns the space of equivariant linear maps from this subrepresentation to another subrepresentation
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``X * self.image(g) = repR.image(g) * X``
        %
        % Example:
%         %   >>> S3 = replab.S(3);
%         %   >>> rep = S3.naturalRep;
%         %   >>> Xrep = randn(3, 3);
%         %   >>> triv = rep.subRep([1;1;1]);
%         %   >>> std = rep.maschke(triv);
%         %   >>> E = triv.subEquivariantTo(std);
%         %   >>> Xsub = E.projectFromParent(Xrep);
%         %   >>> isequal(size(Xsub), [2 1]) % map to std (dim=1) from triv (dim=1)
%         %       1
%         %   >>> norm(Xsub) <= 1e-15
%         %       1
        %
        % Args:
        %   repR (`+replab.SubRep`): Subrepresentation on the target/row space
        %
        % Keyword Args:
        %   special (charstring, optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %   parent (`.Equivariant`, optional): Equivariant space of ``repR.parent`` and ``self.parent``, default: ``[]``
        %
        % Returns:
        %   `+replab.SubEquivariant`: The equivariant vector space
            E = replab.SubEquivariant.make(repR, self, varargin{:});
        end

    end

    methods (Access = protected)

        function c = computeInjectionConditionNumberEstimate(self)
            if all(all(self.injection_internal == self.projection_internal'));
                c = 1;
            else
                c = replab.numerical.condUpperBound(self.injection_internal, self.projection_internal);
            end
        end

        function e = computeProjectorErrorBound(self)
            PiA = self.projector('double/sparse');
            [S DS] = self.parent.commutant.project(PiA);
            e = norm(S - PiA, 'fro'); % + DS; % TODO
        end

    end

    methods (Access = protected) % Implementations

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.SubRep(newParents{1}, self.injection_internal, self.projection_internal);
        end

        function rho = image_double_sparse(self, g)
            if self.isExact && replab.init.cyclolab().works
                rho = double(self.image_exact(g));
            else
                rho = self.projection('double/sparse') * self.parent.image(g, 'double/sparse') * self.injection('double/sparse');
            end
        end

        function rho = image_exact(self, g)
            assert(self.isExact);
            rho = self.projection('exact') * self.parent.image(g, 'exact') * self.injection('exact');
        end

        function c = computeConditionNumberEstimate(self)
            if self.isUnitary
                c = 1;
            else
                % rho = parent
                % rhoU is the unitarization of rho with rho = B * rhoU * Binv
                % self = P B rhoU Binv I
                c = self.injectionConditionNumberEstimate * self.parent.conditionNumberEstimate;
            end
        end

        function e = computeErrorBound(self)
            if self.isSimilarRep
                I = self.injection('double/sparse');
                P = self.projection('double/sparse');
                d = self.dimension;
                prodError = self.biorthogonalityErrorBound; % || dP I ||F
                % let P = I^-1
                % we assume I is exact, and P~ = projection_internal, with P~ = P + dP and dP the error
                % || P rho I - P~ rho~ I ||F =~ || P drho I ||F + || P~ rho I ||F =
                % = ||P||2 ||drho||F ||I||2 + || dP I P rho I||F
                %
                % now || dP I P rho I||F = || dP I ||F ||P rho I||2
                %
                % thus e = cond(injection) * parent.errorBound + self.conditionNumberEstimate*prodError
                term1 = self.injectionConditionNumberEstimate * self.parent.errorBound; % ||P||2 ||I||2 * ||drho||F
                term2 = self.conditionNumberEstimate * prodError;
                e = term1 + term2;
            else
                eU = self.parent.conditionNumberEstimate * self.projectorErrorBound;
                dU = pi*sqrt(self.dimension/2)*eU/(1-eU);
                e1 = self.injectionConditionNumberEstimate * self.parent.errorBound; % parent error
                e2 = self.injectionConditionNumberEstimate * self.parent.conditionNumberEstimate * (2*dU + dU^2);
                e = e1 + e2;
            end
        end

        %function [A Ainv] = unitaryChangeOfBasis(self) TODO restore
        %    if self.parent.knownUnitary
        %        P = self.projection('double/sparse');
        %        X = P * P';
        %        A = chol(X, 'lower');
        %        Ainv = inv(A);
        %    else
        %        [A Ainv] = unitaryChangeOfBasis@replab.Rep(self);
        %    end
        %end

    end
% $$$         function verifyInvariance(self)
% $$$         % Verifies that the basis of this object defines a subrepresentation up to a given precision
% $$$         %
% $$$         % From the basis `.basis` which we write $B$ of this subrepresentation, we compute a projector
% $$$         % $\pi_R = B B^\dagger$. The user provides a parameter ``epsilon`` ($\epsilon$), and when the verification
% $$$         % succeeds, it means that there exists a projector $\pi_K$ on an invariant subspace
% $$$         % with $|| \pi_R - \pi_K ||_2 \le \epsilon$.
% $$$         %
% $$$         % Limitations:
% $$$         %
% $$$         % * Assumes that the parent representation is unitary, and the basis is
% $$$         %   orthonormal.
% $$$         % * Assumes that the errors in the basis are much larger than the errors
% $$$         %   in the representation images themselves.
% $$$         % * Only works for finite groups.
% $$$         %
% $$$         % Args:
% $$$         %   epsilon (double): Maximal error on the subrepresentation basis (see above)
% $$$             assert(isa(self.group, 'replab.FiniteGroup'));
% $$$             assert(isequal( [isUnitary] )); % ETC
% $$$             % by the construction of the BSGS chain, there is a family of sets {T_i},
% $$$             % where i = 1,...,depth, such that every group element can be written as
% $$$             % g = t_1 t_2 ... t_depth, and t_i \in T_i
% $$$             %
% $$$             % As S, we take the union of the T_i. Done.
% $$$             %
% $$$             % Then k = depth, delta = 0, q = 0
% $$$             %
% $$$             % Use Algorithm IV.1
% $$$         end

% $$$         function verifyIrreducibilityExact(self)
% $$$         % Verifies that this subrepresentation is irreducible
% $$$         %
% $$$         % * Assumes that the basis is exact
% $$$         % * Same assumptions as `.verifyInvariance`
% $$$         %
% $$$         % Args:
% $$$         %   pthr (double): Probability of a false positive (should be ~1e-12 or something)
% $$$         %   thetamax (double): Upper bound on theta
% $$$             assert(isa(self.group, 'replab.FiniteGroup'));
% $$$             % S = the whole group, and we sample from the Haar measure
% $$$             % t = depth = 1, choose m such that the theta is << 1
% $$$             % quality of the generating set enters in the length of the product
% $$$             %
% $$$             % Estimate the probability of a false negative given the above
% $$$             % (run experimental experiments?)
% $$$         end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'injection_internal';
            names{1, end+1} = 'projection_internal';
            names{1, end+1} = 'mapsAreAdjoint';
            names{1, end+1} = 'isSimilarRep';
        end

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.Rep(self);
            if self.dimension < 15
                if self.injectionIsExact && replab.init.cyclolab().works
                    type = 'exact';
                else
                    type = 'double';
                end
                for i = 1:self.dimension
                    names{1, end+1} = sprintf('basis(%d,''%s'')', i, type);
                    values{1, end+1} = full(self.basis(i, type));
                end
            end
        end

        % Obj

        function l = laws(self)
            l = replab.laws.SubRepLaws(self);
        end

        % Rep

        function b = isExact(self)
            b = self.mapsAreExact && self.parent.isExact;
        end

        function rep = complexification(self)
            rep = self.parent.complexification.subRep(self.injection_internal, 'projection', self.projection_internal);
        end

        function rep = dual(self)
            if self.knownIrreducible
                args = {'isIrreducible', true};
            else
                args = {};
            end
            rep = self.parent.dual.subRep(self.projection_internal.', 'projection', self.injection_internal.', args{:});
        end

        function rep = conj(self)
            if self.knownIrreducible
                args = {'isIrreducible', true};
            else
                args = {};
            end
            rep = self.parent.conj.subRep(conj(self.injection_internal), 'projection', conj(self.projection_internal), args{:});
        end

        % Rep: Equivariant spaces

        function a = antilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            a = self.cached(['antilinearInvariant_' type], @() self.subEquivariantFrom(conj(self),  'special', 'antilinear', 'type', type));
        end

        function b = bilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            b = self.cached(['bilinearInvariant_' type], @() self.subEquivariantTo(dual(self),  'special', 'bilinear', 'type', type));
        end

        function c = commutant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            c = self.cached(['commutant_' type], @() self.subEquivariantFrom(self,  'special', 'commutant', 'type', type));
        end

        function h = hermitianInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            h = self.cached(['hermitianInvariant_' type], @() self.subEquivariantTo(conj(dual(self)),  'special', 'hermitian', 'type', type));
        end

        function h = sesquilinearInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            h = self.cached(['sesquilinearInvariant_' type], @() self.subEquivariantTo(conj(dual(self)),  'special', 'sesquilinear', 'type', type));
        end

        function h = symmetricInvariant(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            h = self.cached(['symmetricInvariant_' type], @() self.subEquivariantTo(dual(self),  'special', 'symmetric', 'type', type));
        end

        function t = trivialRowSpace(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            t = self.cached(['trivialRowSpace_' type], @() self.subEquivariantTo(...
                replab.SubRep.identical(self.group.trivialRep(self.field, self.dimension)), ...
                'special', 'trivialRows', 'type', type));
        end

        function t = trivialColSpace(self, type)
            if nargin < 2 || isempty(type) || strcmp(type, 'double/sparse')
                type = 'double';
            end
            t = self.cached(['trivialColSpace_' type], @() self.equivariantFrom(...
                replab.SubRep.identical(self.group.trivialRep(self.field, self.dimension)), ...
                'special', 'trivialCols', 'type', type));
        end

    end

    methods (Static)

        function sub = identical(parent)
        % Creates a full subrepresentation of the given representation, with identity projection/injection maps
        %
        % Args:
        %   parent (`+replab.Rep`): Representaiton
        %
        % Returns:
        %   `+replab.SubRep`: Subrepresentation identical to ``parent``
            d = parent.dimension;
            sub = parent.subRep(speye(d), 'projection', speye(d));
            sub.copyProperties(parent);
        end

        function subRep = directSumBiorthogonal(parent, subReps)
        % Computes the direct sum of subrepresentations of the same parent representation
        %
        % The subrepresentations must have biorthogonal injection projection maps:
        % ``subReps{i}.projection * subReps{i}.injection ~= eye(subReps{i}.dimension)``, but
        % ``subReps{i}.projection * subReps{j}.injection ~= 0`` for ``i ~= j``.
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation
        %   subReps (cell(1,\*) of `+replab.SubRep`): Subrepresentations of the parent representation
        %
        % Returns:
        %   `+replab.SubRep`: A block-diagonal subrepresentation composed of the given subrepresentations
            Ps = cellfun(@(sr) sr.projection_internal, subReps, 'uniform', 0);
            Is = cellfun(@(sr) sr.injection_internal, subReps, 'uniform', 0);
            newI_internal = horzcat(Is{:});
            newP_internal = vertcat(Ps{:});
            subRep = parent.subRep(newI_internal, 'projection', newP_internal);
        end

    end

end
