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
        hasExactBasis % (logical): Whether the change of basis matrices are exact
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
            A_unitary = all(all(A_internal' == Ainv_internal));
            if A_unitary && parent.isUnitary
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'isUnitary', true);
                assert(~exists || isempty(oldValue) || isequal(oldValue, true), 'This representation is actually unitary');
            end
            if parent.inCache('trivialDimension')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'trivialDimension', parent.trivialDimension);
                assert(~exists || isempty(oldValue) || oldValue ~= parent.trivialDimension);
            end
            if parent.inCache('isIrreducible')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'isIrreducible', parent.isIrreducible);
                assert(~exists || isempty(oldValue) || oldValue ~= parent.isIrreducible);
            end
            if parent.inCache('frobeniusSchurIndicator')
                [restArgs, exists, oldValue] = replab.util.keyValuePairsUpdate(restArgs, 'frobeniusSchurIndicator', parent.frobeniusSchurIndicator);
                assert(~exists || isempty(oldValue) || oldValue ~= parent.frobeniusSchurIndicator);
            end
            self@replab.Rep(parent.group, parent.field, d, restArgs{:});
            % populate SimilarRep properties
            if ~parent.isExact
                A_internal = double(A_internal);
                Ainv_internal = double(Ainv_internal);
            end
            if isa(A_internal, 'replab.cyclotomic') && isa(Ainv_internal, 'replab.cyclotomic')
                prodError = 0;
            else
                prodError = double(norm(A_internal*Ainv_internal - eye(d), 'fro'));
            end
            isExact = prodError == 0;
            if ~isExact
                A_internal = double(A_internal);
                Ainv_internal = double(Ainv_internal);
            end
            self.parent = parent;
            self.A_internal = A_internal;
            self.Ainv_internal = Ainv_internal;
            self.A_Ainv_error = prodError;
            if isempty(basisConditionNumberEstimate)
                if A_unitary
                    basisConditionNumberEstimate = 1;
                else
                    basisConditionNumberEstimate = replab.numerical.condUpperBound(A_internal, Ainv_internal);
                end
            end
            self.basisConditionNumberEstimate = basisConditionNumberEstimate;
            self.hasExactBasis = isExact;
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

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = 'Similar representation';
        end

        % Rep

        function b = isExact(self)
            b = self.hasExactBasis && self.parent.isExact;
        end

    end


    methods (Access = protected) % Implementations

        % Rep

        function rho = image_double_sparse(self, g)
            if self.isExact
                rho = double(self.image_exact(g));
            else
                rho = self.A_internal * self.parent.image(g, 'double/sparse') * self.Ainv_internal;
            end
        end

        function rho = image_exact(self, g)
            rho = self.A_internal * self.parent.image(g, 'exact') * self.Ainv_internal;
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
            if self.cachedOrDefault('isUnitary', false)
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
            res = replab.SimilarRep(rep, speye(d), speye(d), 'isUnitary', rep.isUnitary);
            if rep.inCache('isDivisionAlgebraCanonical')
                res.cache('isDivisionAlgebraCanonical', rep.isDivisionAlgebraCanonical, 'error');
            end
        end

    end

end
