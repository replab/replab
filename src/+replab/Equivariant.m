classdef Equivariant < replab.Domain
% Describes a vector space of group-equivariant matrices
%
% Let ``repC`` and ``repR`` be two representations of the same group ``G``.
%
% This describes the set of matrices ``X`` such that ``repR.image(g) * X = X * repC.image(g)``
%
% See Proposition 4 of
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
%
% There are two special cases of equivariant spaces.
%
% - ``commutant`` equivariant spaces have ``repC == repR``,
% - ``hermitian`` equivariant spaces have ``repC == repR.conjugate.dual``.
% - ``trivialRows`` equivariant spaces have ``repR == repC.group.trivialRep(field, repR.dimension)``
% - ``trivialCols`` equivariant spaces have ``repC == repR.group.trivialRep(field, repC.dimension)``
%
% When ``repR`` is unitary, the ``commutant`` and ``hermitian`` cases are identical.

    properties (SetAccess = protected)
        field % ({'R', 'C'}): Field of the vector space real (R) or complex x(C)
        nR % (integer): Row size
        nC % (integer): Column size
        group % (`+replab.CompactGroup`): Group being represented
        repR % (`+replab.Rep`): Representation of row space
        repC % (`+replab.Rep`): Representation of column space
        special % ({'hermitian', 'commutant', 'trivialRows', 'trivialCols', []}): Whether the equivariant space has special structure
    end

    properties (Access = protected)
        domain % (`+replab.Domain`): Domain, real or complex matrices
        cachedSamples_ % (struct of cell(1,\*) of domain elements): Samples
        cachedErrors_ % (struct of double(1,\*)): Error information
    end

    methods (Access = protected)

        function self = Equivariant(repR, repC, special)
        % Constructor; use `+replab.Equivariant.make` in user code
        %
        % That function selects an optimized implementation depending on the use case.
            self.repR = repR;
            self.nR = repR.dimension;
            self.repC = repC;
            self.nC = repC.dimension;
            assert(isequal(repR.field, repC.field), ...
                   'Both representations must have be defined on the same field');
            assert(isempty(special) || ismember(special, {'hermitian', 'commutant', 'trivialRows', 'trivialCols'}));
            self.field = repR.field;
            assert(repR.group == repC.group, ...
                   'Both representations must be defined on the same group');
            self.group = repR.group;
            self.domain = replab.domain.Matrices(self.field, self.nR, self.nC);
            self.special = special;
            self.cachedSamples_ = struct;
            self.cachedErrors_ = struct;
        end

    end

    methods (Access = protected) % Protected methods

        function X1 = project_exact(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Implementation of `.project`
        %
        % Raises:
        %   An error if `.isExact` is false.
        %
        % Args:
        %   X (`.cyclotomic`(\*,\*)): Matrix to project
        %
        % Returns
        % -------
        %   X1: `.cyclotomic`(\*,\*)
        %     Projected matrix
            error('Exact projection not implemented');
        end

        function X1 = project_intval(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Implementation of `.project`
        %
        % Raises:
        %   An error if `.isExact` is false.
        %
        % Args:
        %   X (intval(\*,\*)): Matrix to project
        %
        % Returns
        % -------
        %   X1: `.intval`(\*,\*)
        %     Projected matrix
            error('Projection in interval arithmetic not implemented');
        end

        function [X1, err] = project_double_sparse(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Implementation of `.project`
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix to project
        %
        % Returns
        % -------
        %   X1: double(\*,\*)
        %     Projected matrix
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned ``X1`` to the invariant subspace in Frobenius norm
            error('Abstract');
        end

    end

    methods % Projection

        function [X1, err] = project(self, X, type)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Performs the integration
        %
        % `` X1 = int{g in G} dg rhoR.image(g) * X * rhoC.inverseImage(g) ``
        %
        % Note that there are little benefits usually to use the ``'double/sparse'`` type compared to ``'double'``.
        %
        % Raises:
        %   An error if ``type`` is ``'exact'`` and `.isExact` is false.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   X1: double(\*,\*) or `.cyclotomic`(\*,\*)
        %     Projected matrix
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned ``X1`` to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'double'
                if nargout > 1
                    [X1, err] = self.project_double_sparse(X);
                else
                    X1 = self.project_double_sparse(X);
                end
                X1 = full(X1);
              case 'double/sparse'
                if nargout > 1
                    [X1, err] = self.project_double_sparse(X);
                else
                    X1 = self.project_double_sparse(X);
                end
              case 'intval'
                if ~isa(X, 'intval')
                    X = intval(X);
                end
                X1 = self.project_intval(X);
              case 'exact'
                if isa(X, 'double')
                    X = replab.cyclotomic.fromDoubles(X);
                end
                X1 = self.project_exact(X);
              otherwise
                error('Type must be either double, double/sparse or exact');
            end
        end

    end

    methods

        function b = isExact(self)
        % Returns whether this equivariant space can compute exact projection
        %
        % Returns:
        %   logical: True if the call ``self.projection(X, 'exact')`` always succeeds
            b = false;
        end

        function b = hasIntval(self)
        % Returns whether this equivariant space can compute projection in interval arithmetic
        %
        % Returns:
        %   logical: True if the call ``self.projection(X, 'intval')`` always succeeds
            b = false;
        end

        function E1 = subEquivariant(self, subR, subC, varargin)
        % Constructs a invariant subspace of an equivariant space
        %
        % Args:
        %   subC (`+replab.SubRep`): A subrepresentation of ``self.repC``
        %   subR (`+replab.SubRep`): A subrepresentation of ``self.repR``
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
            assert(isa(subC, 'replab.SubRep'));
            assert(isa(subR, 'replab.SubRep'));
            assert(subC.parent == self.repC);
            assert(subR.parent == self.repR);
            E = subR.subEquivariantFrom(subC, 'parent', self, varargin{:});
        end

    end

    methods % Cached samples

        function clearCache(self, context)
        % Clears the samples cached for the given context
        %
        % Args:
        %   context (`+replab.Context`): Context to clear
            self.cachedSamples_ = rmfield(self.cachedSamples_, context.id);
            self.cachedErrors_ = rmfield(self.cachedErrors_, context.id);
        end

        function [X err] = sampleInContext(self, context, ind)
        % Returns an approximate sample from this equivariant space along with estimated numerical error
        %
        % The samples are cached in a context.
        %
        % Args:
        %   context (`+replab.Context`): Context in which samples are cached
        %   ind (double): 1-based index of the sample
        %
        % Returns
        % -------
        % X:
        %   double(\*,\*): A sample from this equivariant space
        % err:
        %   double: Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %           the invariant subspace in Frobenius norm
            assert(~context.closed);
            id = context.id;
            if ~isfield(self.cachedSamples_, id)
                context.register(self);
                self.cachedSamples_.(id) = cell(1, 0);
                self.cachedErrors_.(id) = zeros(1, 0);
            end
            n = length(self.cachedSamples_.(id));
            samples = self.cachedSamples_.(id);
            errors = self.cachedErrors_.(id);
            if ind > n
                for i = n+1:ind
                    [X err] = self.sample;
                    samples{1, i} = X;
                    errors{1, i} = err;
                end
                self.cachedSamples_.(id) = samples;
                self.cachedErrors_.(id) = errors;
            end
            X = samples{ind};
            err = errors(ind);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            if isequal(self.special, 'hermitian')
                s = sprintf('%d x %d invariant Hermitian matrices over %s', ...
                            self.nR, self.nC, self.field);
            elseif isequal(self.special, 'commutant')
                s = sprintf('%d x %d commutant matrices over %s', ...
                            self.nR, self.nC, self.field);
            else
                s = sprintf('%d x %d equivariant matrices over %s', ...
                            self.nR, self.nC, self.field);
            end
        end

        % Domain

        function l = laws(self)
            l = replab.laws.EquivariantLaws(self);
        end

        function b = eqv(self, X, Y)
            b = self.domain.eqv(X, Y);
        end

        function [X err] = sample(self, type)
        % Returns an approximate sample from this equivariant space along with estimated numerical error
        %
        % Raises:
        %   An error if ``type`` is ``'exact'`` and `.isExact` is false.
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   X: double(\*,\*)
        %     A sample from this equivariant space
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned ``X`` to
        %     the invariant subspace in Frobenius norm
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            if strcmp(type, 'exact')
                i = randi(self.nR);
                j = randi(self.nC);
                X = replab.cyclotomic.sparse(i, j, 1, self.nR, self.nC);
                X = self.project(X, 'exact');
                err = 0;
            else
                [X err] = self.project(self.domain.sample, type);
            end
        end

    end

    methods (Static)

        function E = make(repR, repC, varargin)
        % Returns the space of equivariant linear maps between two representations
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`+replab.Rep`): Representation on the target/row space
        %   repC (`+replab.Rep`): Representation on the source/column space
        %
        % Keyword Args:
        %   special ('commutant', 'hermitian', 'trivialRows', 'trivialCols' or '', optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `+replab.Equivariant`: The equivariant vector space
            args = struct('special', '', 'type', 'double');
            args = replab.util.populateStruct(args, varargin);
            if isa(repR, 'replab.SimilarRep') && isa(repC, 'replab.SimilarRep')
                E = replab.equi.Equivariant_forSimilarRep.make(repR, repC, varargin{:});
            elseif isa(repR, 'replab.SimilarRep') && isa(repC, 'replab.SubRep')
                E = replab.Equivariant.make(repR.toSubRep, repC, varargin{:});
            elseif isa(repR, 'replab.SubRep') && isa(repC, 'replab.SimilarRep')
                E = replab.Equivariant.make(repR, repC.toSubRep, varargin{:});
            elseif isa(repR, 'replab.SubRep') && isa(repC, 'replab.SubRep')
                E = replab.SubEquivariant.make(repR, repC, varargin{:});
            elseif isa(repR.group, 'replab.FiniteGroup')
                if repR.group.order < 65536 && (~repR.isExact || ~repC.isExact)
                    E = replab.equi.Equivariant_forFiniteGroup_explicitSum(repR, repC, args.special);
                else
                    E = replab.equi.Equivariant_forFiniteGroup_relativeReynolds(repR, repC, args.special);
                end
            else
                E = replab.equi.Equivariant_forCompactGroup(repR, repC, args.special);
            end
        end

    end

end
