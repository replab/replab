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

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation of dimension $D$
        injection_internal % (double(D,d) or `.cyclotomic`(D,d), may be sparse): Injection map
        projection_internal % (double(d,D) or `.cyclotomic`(d,D), may be sparse): Projection map
        hasExactMaps % (logical): Whether the injection and projection maps are exact
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
            assert(size(projection_internal, 1) == d);
            assert(size(projection_internal, 2) == D);
            assert(parent.dimension == D, 'Incorrect dimension');
            if parent.overR
                assert(isreal(injection_internal) && isreal(projection_internal), 'A real Rep can only have real subrepresentations.');
            end
            args = struct('injectionConditionNumberEstimate', [], 'projectorErrorBound', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            IP_unitary = all(all(injection_internal == projection_internal'));
            if IP_unitary && parent.knownUnitary
                restArgs = replab.util.keyValuePairsUpdate(restArgs, 'knownUnitary', true);
            end
            if parent.inCache('trivialDimension') && parent.trivialDimension == 0
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', 0);
                assert(~exists || oldValue == 0);
            end
            self@replab.Rep(parent.group, parent.field, d, restArgs{:});
            if isa(injection_internal, 'replab.cyclotomic') && isa(projection_internal, 'replab.cyclotomic')
                hasExactMaps = true;
                assert(isempty(args.projectorErrorBound) || args.projectorErrorBound == 0);
                args.projectorErrorBound = 0;
            elseif replab.numerical.isExact(injection_internal) && replab.numerical.isExact(projection_internal) && ...
                    ~isempty(args.projectorErrorBound) && args.projectorErrorBound == 0
                hasExactMaps = true;
            else
                hasExactMaps = false;
            end
            self.parent = parent;
            self.injection_internal = injection_internal;
            self.projection_internal = projection_internal;
            self.hasExactMaps = hasExactMaps;
            if ~isempty(args.projectorErrorBound)
                self.cache('projectorErrorBound', args.projectorErrorBound, 'error');
            end
            if ~isempty(args.injectionConditionNumberEstimate)
                self.cache('injectionConditionNumberEstimate', args.injectionConditionNumberEstimate, 'error');
            end
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_SubRepOfSubRep(self, options)
            if isa(self.parent, 'replab.SubRep')
                if self.isIntegerValued || self.parent.isIntegerValued || ...
                        (options.dense && (options.approximate || (self.hasExactMaps && self.parent.hasExactMaps)))
                    newI = self.parent.injection_internal * self.injection_internal;
                    newP = self.projection_internal & self.parent.projection_internal;
                    res = replab.SubRep(self.parent.parent, newI, newP);
                    return
                end
            end
            res = [];
        end

        function res = rewriteTerm_SubRepOfSimilarRep(self, options)
            if isa(self.parent, 'replab.SimilarRep')
                if self.isIntegerValued || self.parent.isIntegerValued || ...
                        (options.dense && (options.approximate || (self.hasExactMaps && self.parent.hasExactBasis)))
                    newI = self.parent.Ainv_internal * self.injection_internal;
                    newP = self.projection_internal * self.parent.A_internal;
                    res = replab.SubRep(self.parent.parent, newI, newP);
                    return
                end
            end
            res = [];
        end

    end

    methods

        function sub1 = withNoise(self, injectionMapNoise, projectionMapNoise)
        % Adds Gaussian noise to the injection/projection maps
        %
        % This method has two call conventions.
        %
        % If ``sub.withNoise(sigma)`` is called, then the injection map is changed to
        % ``injection1 = injection + Delta``, where ``Delta`` is a matrix with normally distributed
        % entries of standard deviation ``sigma``. ``Delta`` is real (resp. complex) if the representation
        % is real (resp. complex). Then the original `.projection` map is ignored and ``projection1`` is recovered
        % as in `.Rep.subRep`.
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
                projection1 = inv(projection1 * injection1) * projection1;
                sub1 = self.parent.subRep(injection1, 'projection', projection1);
              otherwise
                error('Wrong calling convention');
            end
        end

        function b = isIntegerValued(self)
        % Returns whether both the injection map and the projection map are expressed with Gaussian integer coefficients
        %
        % Returns:
        %   logical: True if both `.injection` and `.projection` have Gaussian integer entries
            I = self.injection_internal;
            P = self.projection_internal;
            b = self.hasExactMaps && all(all(I == round(I))) && all(all(P == round(P)));
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
            if isempty(indices)
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

        function sub1 = refine(self, varargin)
        % Refines this subrepresentation
        %
        % Assumes that the subrepresentation is already close to an exact subrepresentation, and refines its subspace
        % by an iterative procedure applied on its `.injection` and `.projection` maps.
        %
        % Keyword Args:
        %   largeScale (logical or ``[]``, optional): Whether to use the large-scale version of the algorithm, default ``[]`` (automatic selection)
        %   numNonImproving (integer, optional): Number of non-improving steps before stopping the large-scale algorithm, default ``3``
        %   nSamples (integer, optional): Number of samples to use in the large-scale version of the algorithm, default ``3``
        %   nInnerIterations (integer, optional): Number of inner iterations in the medium-scale version of the algorithm, default ``10``
        %   maxIterations (integer, optional): Maximum number of (outer) iterations, default ``1000``
        %
        % Returns
        % -------
        %   sub1: `.SubRep`
        %     Subrepresentation with refined subspace (injection/projection maps)
            args = struct('numNonImproving', 5, 'largeScale', self.parent.dimension > 1000, 'nSamples', 3, 'nInnerIterations', 10, 'maxIterations', 1000);
            args = replab.util.populateStruct(args, varargin);
            I = self.injection;
            P = self.projection;
            if self.parent.knownUnitary
                if ~all(all(self.injection == self.projection'))
                    [I, ~] = qr(I, 0);
                end
                if args.largeScale
                    I1 = replab.rep.refine_unitaryLargeScale(self.parent, I, args.numNonImproving, args.nSamples, args.maxIterations);
                else
                    I1 = replab.rep.refine_unitaryMediumScale(self.parent, I, args.nInnerIterations, args.maxIterations);
                end
                sub1 = self.parent.subRep(I1, 'projection', I1', 'isUnitary', true);
            else
                if args.largeScale
                    [I1, P1] = replab.rep.refine_nonUnitaryLargeScale(self.parent, I, P, args.numNonImproving, args.nSamples, args.maxIterations);
                else
                    [I1, P1] = replab.rep.refine_nonUnitaryMediumScale(self.parent, I, P, args.nInnerIterations, args.maxIterations);
                end
                sub1 = self.parent.subRep(I1, 'projection', P1);
            end
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
            e = norm(S - PiA, 'fro'); % + DS;
        end

    end

    methods (Access = protected) % Implementations

        function r = computeDouble(self)
            args = self.knownProperties({'injectionConditionNumberEstimate'});
            r = replab.SubRep(double(self.parent), self.injection('double/sparse'), self.projection('double/sparse'), args{:});
        end

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.SubRep(newParents{1}, self.injection_internal, self.projection_internal);
        end

        function rho = image_double_sparse(self, g)
            if self.isExact
                rho = double(self.image_exact(g));
            else
                rho = self.projection_internal * self.parent.image(g, 'double/sparse') * self.injection_internal;
            end
        end

        function rho = image_exact(self, g)
            assert(self.isExact);
            rho = self.projection_internal * self.parent.image(g, 'exact') * self.injection_internal;
        end

        function c = computeConditionNumberEstimate(self)
            if self.cachedOrDefault('isUnitary', false)
                c = 1;
            else
                % rho = parent
                % rhoU is the unitarization of rho with rho = B * rhoU * Binv
                % self = P B rhoU Binv I
                c = self.injectionConditionNumberEstimate * self.parent.conditionNumberEstimate;
            end
        end

        function e = computeErrorBound(self)
            eU = self.parent.conditionNumberEstimate * self.projectorErrorBound;
            dU = pi*sqrt(self.dimension/2)*eU/(1-eU);
            e1 = self.injectionConditionNumberEstimate * self.parent.errorBound; % parent error
            e2 = self.injectionConditionNumberEstimate * self.parent.conditionNumberEstimate * (2*dU + dU^2);
            e = e1 + e2;
        end

        function [A Ainv] = unitaryChangeOfBasis(self)
            if self.parent.knownUnitary
                P = self.projection('double/sparse');
                X = P * P';
                A = chol(X, 'lower');
                Ainv = inv(A);
            else
                [A Ainv] = unitaryChangeOfBasis@replab.Rep(self);
            end
        end

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

% $$$         function [s better] = nice(self)
% $$$         % Returns a representation similar to the current subrepresentation, with a nicer basis
% $$$         %
% $$$         % The "niceness" of the basis is implementation dependent. As of the first implementation
% $$$         % of this feature, RepLAB tries to make the basis real, and then with small integer
% $$$         % coefficients.
% $$$         %
% $$$         % The returned subrepresentation is not necessarily unitary.
% $$$         %
% $$$         % In the case no improvement could be made, the original subrepresentation is returned.
% $$$         %
% $$$         % Returns:
% $$$         %   `+replab.SubRep`: A subrepresentation of ``self.parent``
% $$$             s = replab.nice.niceSubRep(self);
% $$$         end


    methods % Implementations

        % Str

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.dimension < 15
                if self.isExact
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
            b = self.hasExactMaps && self.parent.isExact;
        end

    end

    methods (Static)

        function sub = fullSubRep(parent)
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
            subRep = parent.subRep(newI_internal, newP_internal);
        end

    end

end
