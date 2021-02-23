classdef SimilarRep < replab.Rep
% Describes a representation similar to a given representation under a change of basis
%
% We use left action convention, which means that ``image(g) = A * parent.image(g) * inv(A)``,
% where `.A` is the change of basis matrix.
%
% The type of the change of basis matrix can be either double floating-point (and can be sparse),
% or can be given as a `.cyclotomic` matrix; its inverse must be provided as well.
%
% If the following conditions are met:
%
% * the parent representation is exact,
% * the product of the given change of basis matrix and the given inverse is exactly the identity,
%   regardless of its storage as cyclotomics or floating point numbers,
%
% then this representation is exact as well.
%
% Otherwise this representation is inexact.

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Representation this is derived from
        A_internal % (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Change of basis matrix
        Ainv_internal % (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Inverse of change of basis matrix
        A_Ainv_error % (double): Frobenius norm of ``A * Ainv - identity``
        basisConditionNumberEstimate % (double): Estimate of the condition number of `.A_internal`
        mapsAreAdjoint % (logical): True if `.parent` is known to be unitary and `.A_internal` is the conjugate transpose of `.Ainv_internal`
    end

    methods

        function self = SimilarRep(parent, A_internal, Ainv_internal, varargin)
        % Constructs a finite dimension similar to a parent representation
        %
        % Additional keyword arguments are passed to the `.Rep` constructor.
        %
        % If the parent representation is not exact, the change of basis matrices will be converted to
        % floating-point matrices.
        %
        % The `.trivialDimension`, `.isIrreducible`, `.frobeniusSchurIndicator` properties are
        % preserved from the parent. The keyword argument ``isUnitary`` must be provided, except if
        % the parent representation is unitary, and the change of basis matrix is unitary as well.
        %
        % Args:
        %   parent % (`+replab.Rep`): Representation this is derived from
        %   A_internal % (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Change of basis matrix
        %   Ainv_internal % (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Approximate inverse of change of basis matrix
        %
        % Keyword Args:
        %   isUnitary (logical, optional): Whether the resulting representation is unitary; may be omitted (see above)
        %   basisConditionNumberEstimate (double): Upper bound of the condition number of ``A_internal``
            if parent.overR
                assert(isreal(A_internal) && isreal(Ainv_internal), 'A real Rep can only be conjugated by a real orthonormal matrix');
            end
            d = parent.dimension;
            assert(size(A_internal, 1) == d && size(A_internal, 2) == d);
            assert(size(Ainv_internal, 1) == d && size(Ainv_internal, 2) == d);
            % populate Rep properties
            args = struct('basisConditionNumberEstimate', []);
            [args, restArgs] = replab.util.populateStruct(args, varargin);
            basisConditionNumberEstimate = args.basisConditionNumberEstimate;
            mapsAreAdjoint = all(all(A_internal' == Ainv_internal));
            if mapsAreAdjoint && parent.isUnitary
                restArgs = replab.util.keyValuePairsUpdate(restArgs, 'isUnitary', true);
            end
            if parent.inCache('trivialDimension')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', parent.trivialDimension);
                assert(~exists || isempty(oldValue) || oldValue == parent.trivialDimension);
            end
            if parent.inCache('isIrreducible')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'isIrreducible', parent.isIrreducible);
                assert(~exists || isempty(oldValue) || oldValue == parent.isIrreducible);
            end
            if parent.inCache('frobeniusSchurIndicator')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'frobeniusSchurIndicator', parent.frobeniusSchurIndicator);
                assert(~exists || isempty(oldValue) || oldValue == parent.frobeniusSchurIndicator);
            end
            self@replab.Rep(parent.group, parent.field, d, restArgs{:});
            % populate SimilarRep properties
            prodError = norm(double(A_internal)*double(Ainv_internal) - eye(d), 'fro');
            self.parent = parent;
            self.A_internal = A_internal;
            self.Ainv_internal = Ainv_internal;
            self.A_Ainv_error = prodError;
            self.mapsAreAdjoint = mapsAreAdjoint;
            if isempty(basisConditionNumberEstimate)
                if mapsAreAdjoint
                    basisConditionNumberEstimate = 1;
                else
                    basisConditionNumberEstimate = replab.numerical.condUpperBound(A_internal, Ainv_internal);
                end
            end
            self.basisConditionNumberEstimate = basisConditionNumberEstimate;
        end

    end

    methods % Simplification rules

        function res = rewriteTerm_isIdentity(self, options)
            if self.basisIsIntegerValued && all(all(self.A_internal == speye(self.dimension)))
                res = self.parent;
            else
                res = [];
            end
        end

        function res = rewriteTerm_similarRepOfSubRep(self, options)
            if isa(self.parent, 'replab.SubRep')
                if self.basisIsIntegerValued || self.parent.mapsAreIntegerValued || ...
                        (options.dense && (options.approximate || (self.basisIsExact && self.parent.mapsAreExact)))
                    newProjection = self.A_internal * self.parent.projection_internal;
                    newInjection = self.parent.injection_internal * self.Ainv_internal;
                    res = replab.SubRep(self.parent.parent, newInjection, newProjection);
                    return
                end
            end
            res = [];
        end

        function res = rewriteTerm_similarRepOfSimilarRep(self, options)
            if isa(self.parent, 'replab.SimilarRep')
                if self.basisIsIntegerValued || self.parent.basisIsIntegerValued || ...
                        (options.dense && (options.approximate || (self.basisIsExact && self.parent.basisIsExact)))
                    newA = self.A_internal * self.parent.A_internal;
                    newAinv = self.parent.Ainv_internal * self.Ainv_internal;
                    res = replab.SimilarRep(self.parent.parent, newA, newAinv);
                    return
                end
            end
            res = [];
        end

    end

    methods

        function b = basisIsExact(self)
        % Returns whether this similarity transformation is written using a cyclotomic or Gaussian integer basis matrix
        %
        % Returns:
        %   logical: True if both `.A` and `.Ainv` are cyclotomic matrices or have Gaussian integer entries
            X = self.A_internal;
            Y = self.Ainv_internal;
            b = (isa(X, 'cyclotomic') || replab.numerical.isExact(X)) && ...
                (isa(Y, 'cyclotomic') || replab.numerical.isExact(Y));
        end

        function b = basisIsIntegerValued(self)
        % Returns whether this similarity transformation has Gaussian integer change of basis matrices
        %
        % Returns:
        %   logical: True if both `.A` and `.Ainv` have Gaussian integer entries
            A = self.A_internal;
            Ainv = self.Ainv_internal;
            b = self.basisIsExact && all(all(A == round(double(A)))) && all(all(Ainv == round(double(Ainv))));
        end

        function mat = A(self, type)
        % Returns the change of basis matrix
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): Change of basis matrix
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            mat = replab.numerical.convert(self.A_internal, type);
        end

        function mat = Ainv(self, type)
        % Returns the inverse of the change of basis matrix
        %
        % rgs:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or `.cyclotomic`(\*,\*): Inverse of the change of basis matrix
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            mat = replab.numerical.convert(self.Ainv_internal, type);
        end

        function res = collapse(self)
        % Simplifies a SimilarRep of a SubRep/SimilarRep
        %
        % Raises:
        %   An error if `.parent` is not of type `.SubRep` or `.SimilarRep`
        %
        % Returns:
        %   `.Rep`: A subrepresentation of ``.parent.parent`` of the same type as `.parent`
            parent = self.parent;
            switch class(self.parent)
              case 'replab.SubRep'
                newI_internal = parent.injection_internal * self.Ainv_internal;
                newP_internal = self.A_internal * parent.projection_internal;
                res = parent.parent.subRep(newI_internal, 'projection', newP_internal);
                res.copyProperties(self);
              case 'replab.SimilarRep'
                newAinv_internal = parent.Ainv_internal * self.Ainv_internal;
                newA_internal = self.A_internal * parent.A_internal;
                res = parent.parent.similarRep(newA_internal, 'inverse', newAinv_internal);
                res.copyProperties(self);
              otherwise
                error('Not supported');
            end
        end

        function res = toSubRep(self)
        % Returns a subrepresentation of the parent representation equivalent to this representation
        %
        % Returns:
        %   `.SubRep`: A subrepresentation of `.parent`
            res = self.parent.subRep(self.Ainv_internal, 'projection', self.A_internal);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = 'Similar representation';
        end

        % Rep

        function b = isExact(self)
            b = self.basisIsExact && self.parent.isExact;
        end

        function rep = complexification(self)
            rep = self.parent.complexification.similarRep(self.A_internal, 'inverse', self.Ainv_internal);
        end

        function rep = dual(self)
            rep = self.parent.dual.similarRep(self.Ainv_internal.', 'inverse', self.A_internal.');
        end

        function rep = conj(self)
            rep = self.parent.conj.similarRep(conj(self.A_internal), 'inverse', conj(self.Ainv_internal));
        end

    end


    methods (Access = protected) % Implementations

        % Rep

        function c = decomposeTerm(self)
            c = {self.parent};
        end

        function r = composeTerm(self, newParents)
            r = replab.SimilarRep(newParents{1}, self.A_internal, self.Ainv_internal);
        end

        function rho = image_double_sparse(self, g)
            if self.isExact
                rho = double(self.image_exact(g));
            else
                rho = self.A_internal * self.parent.image(g, 'double/sparse') * self.Ainv_internal;
            end
        end

        function rho = image_exact(self, g)
            assert(self.isExact);
            rho = self.A_internal * self.parent.image(g, 'exact') * self.Ainv_internal;
        end

        function b = computeIsIrreducible(self)
            b = self.parent.isIrreducible;
        end

        function e = computeErrorBound(self)

        % let B = A^-1
        % we assume A is exact, and B~ = Ainv_internal, with B~ = B + dB and dB the error
        % || A rho B - A rho~ B~ ||F =~ || A drho B ||F + || A rho dB ||F =
        % = ||A||2 ||drho||F ||B||2 + ||A rho B A dB||F
        %
        % now ||A rho B A dB||F = ||A rho B||2 ||A dB||F
        %
        % thus e = cond(A) * parent.errorBound + self.conditionNumberEstimate*A_Ainv_error
            term1 = self.basisConditionNumberEstimate * self.parent.errorBound; % ||A||2 ||B||2 * ||drho||F
            term2 = self.conditionNumberEstimate * self.A_Ainv_error;
            e = term1 + term2;
        end

        function c = computeConditionNumberEstimate(self)
            if self.isUnitary
                c = 1;
            else
                % rho = parent
                % self = A * rho * Ainv
                %
                % let rhoU be the unitarization of rho, and rho = B * rhoU * Binv
                %
                % then self = A * B * rhoU * Binv * Ainv
                %
                % the condition number of A*B is cond(A) * cond(B)
                c = self.basisConditionNumberEstimate * self.parent.conditionNumberEstimate;
            end
        end

    end

    methods (Static)

        function res = identical(rep)
        % Returns a `.SimilarRep` identical to the given representation
        %
        % Args:
        %   rep (`.Rep`): Representation
        %
        % Returns:
        %   `.SimilarRep`: Identical representation to ``rep``
            d = rep.dimension;
            res = replab.SimilarRep(rep, speye(d), speye(d), 'isUnitary', rep.isUnitary, 'divisionAlgebraName', rep.divisionAlgebraName, args{:});
        end

    end

end
