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
% There are several special cases of equivariant spaces which RepLAB uses extensively; these are parameterized
% by a single representation ``rho`` from which the equivariant space is constructed.
% They are distinguished by the `.special` property.
%
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | Special name   | repR                        | repC                        | Additional constraint   |
% +================+=============================+=============================+=========================+
% | antilinear     | rho                         | conj(rho)                   |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | bilinear       | dual(rho)                   | rho                         |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | commutant      | rho                         | rho                         |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | hermitian      | conj(dual(rho))             | rho                         | X = X'                  |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | sesquilinear   | conj(dual(rho))             | rho                         |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | symmetric      | dual(rho)                   | rho                         | X = X.'                 |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | trivialRows    | trivial (d = rho.dimension) | rho                         |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
% | trivialCols    | rho                         | trivial (d = rho.dimension) |                         |
% +----------------+-----------------------------+-----------------------------+-------------------------+
%
% When ``rho`` is unitary:
%
% * ``commutant`` and ``sesquilinear`` spaces are identical,
% * ``antilinear`` and ``bilinear`` space are identical.
%
% When ``rho`` is real:
%
% * ``antilinear`` and ``commutant`` spaces are identical,
% * ``bilinear`` and ``sesquilinear`` spaces are identical,
% * ``hermitian`` and ``symmetric`` spaces are identical.

    properties (SetAccess = protected)
        field % ({'R', 'C'}): Field of the vector space real (R) or complex x(C)
        nR % (integer): Row size
        nC % (integer): Column size
        group % (`+replab.CompactGroup`): Group being represented
        repR % (`+replab.Rep`): Representation of row space
        repC % (`+replab.Rep`): Representation of column space
        special % (charstring): Whether the equivariant space has special structure, see comment in `+replab.Equivariant`
    end

    properties (Access = protected)
        domain % (`+replab.Domain`): Domain, real or complex matrices
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
            self.field = repR.field;
            self.group = repR.group;
            self.domain = replab.domain.Matrices(self.field, self.nR, self.nC);
            self.special = special;
        end

    end

    methods (Access = protected) % Protected methods

        function D = computeDecomposition(self)
            switch self.special
              case 'antilinear'
                D = self.repR.decomposition.antilinearInvariant;
              case 'bilinear'
                D = self.repC.decomposition.bilinearInvariant;
              case 'commutant'
                D = self.repR.decomposition.commutant;
              case 'hermitian'
                D = self.repC.decomposition.hermitianInvariant;
              case 'sesquilinear'
                D = self.repC.decomposition.sesquilinearInvariant;
              case 'symmetric'
                D = self.repC.decomposition.symmetricInvariant;
              case 'trivialRows'
                D = self.repC.decomposition.trivialRowSpace;
              case 'trivialCols'
                D = self.repR.decomposition.trivialColSpace;
              case ''
                D = self.repR.decomposition.irreducibleEquivariantFrom(self.repC.decomposition);
            end
        end

        function X1 = project_exact(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Implementation of `.project`
        %
        % Raises:
        %   An error if `.isExact` is false.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic` (\*,\*), may be sparse): Matrix to project; if double, should be converted to `.cyclotomic`
        %
        % Returns
        % -------
        %   X1: `.cyclotomic`(\*,\*)
        %     Projected matrix
            error('Exact projection not implemented');
        end

        function [X1, err] = project_double_sparse(self, X)
        % Projects any ``nR x nC`` matrix in the equivariant subspace
        %
        % Implementation of `.project`
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic` (\*,\*), may be sparse): Matrix to project
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
        %   X (double(\*,\*) or `.cyclotomic` (\*,\*), may be sparse): Matrix to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   X1: double(\*,\*) or `.cyclotomic` (\*,\*)
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
              case 'exact'
                X1 = self.project_exact(X);
              otherwise
                error('Type must be either double, double/sparse or exact');
            end
        end

    end

    methods % Properties

        function b = isExact(self)
        % Returns whether this equivariant space can compute exact projection
        %
        % Returns:
        %   logical: True if the call ``self.projection(X, 'exact')`` always succeeds
            b = false;
        end

    end

    methods % Subspaces

        function E1 = subEquivariant(self, subR, subC, varargin)
        % Constructs a invariant subspace of an equivariant space
        %
        % Args:
        %   subC (`+replab.SubRep`): A subrepresentation of ``self.repC``
        %   subR (`+replab.SubRep`): A subrepresentation of ``self.repR``
        %
        % Keyword Args:
        %   special (charstring, optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
            assert(isa(subC, 'replab.SubRep'));
            assert(isa(subR, 'replab.SubRep'));
            assert(subC.parent == self.repC);
            assert(subR.parent == self.repR);
            E = subR.subEquivariantFrom(subC, 'parent', self, varargin{:});
        end

        function D = decomposition(self)
        % Returns the decomposition of this equivariant space
        %
        % Returns:
        %   `.IrreducibleEquivariant`: The equivariant space
            D = self.cached('decomposition', @() self.computeDecomposition);
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            switch self.special
              case 'antilinear'
                  s = sprintf('%d x %d matrices representing an equivariant antilinear map over %s', ...
                              self.nR, self.nC, self.field);
              case 'bilinear'
                  s = sprintf('%d x %d matrices representing an equivariant bilinear form over %s', ...
                              self.nR, self.nC, self.field);
              case 'commutant'
                s = sprintf('%d x %d commutant matrices over %s', self.nR, self.nC, self.field);
              case 'hermitian'
                  s = sprintf('%d x %d matrices representing an equivariant Hermitian form over %s', ...
                              self.nR, self.nC, self.field);
              case 'sesquilinear'
                  s = sprintf('%d x %d matrices representing an equivariant sesquilinear form over %s', ...
                              self.nR, self.nC, self.field);
              case 'symmetric'
                  s = sprintf('%d x %d matrices representing an equivariant symmetric form over %s', ...
                              self.nR, self.nC, self.field);
              otherwise
                s = sprintf('%d x %d equivariant matrices over %s', self.nR, self.nC, self.field);
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

        function [repR, repC] = validateArguments(repR, repC, special)
        % Validate and complete arguments used to define an equivariant space
        %
        % When the ``special`` argument is specified, one can omit either ``repR`` or ``repC`` and it will be completed
        % automatically.
        %
        % Args:
        %   repR (`.Rep` or ``[]``): Representation acting on rows
        %   repC (`.Rep` or ``[]``): Representation acting on columns
        %   special (charstring or ''): One of the equivariant types, see `.Equivariant`
        %
        % Returns
        % -------
        %   repR: `.Rep`
        %     Row representation
        %   repC: `.Rep`
        %     Column representation
            assert(repR.group == repC.group, 'Both representations must be defined on the same group');
            assert(isequal(repR.field, repC.field), 'Both representations must have be defined on the same field');
            if isempty(special)
                special = '';
            end
            valid = {'', 'antilinear', 'commutant', 'hermitian', 'sesquilinear', 'symmetric', 'trivialRows', 'trivialCols'};
            assert(ismember(special, valid));
            if isempty(special)
                assert(~isempty(repR) && ~isempty(repC), 'The row or column representation can be omitted only for structured (special) spaces.');
            else
                assert(~isempty(repR) || ~isempty(repC), 'Only one of the row or column representation can be omitted');
                if isempty(repR)
                    switch special
                      case 'antilinear'
                        repR = conj(repC);
                      case {'bilinear', 'symmetric'}
                        repR = dual(repC);
                      case 'commutant'
                        repR = repC;
                      case {'sesquilinear', 'hermitian'}
                        repR = dual(conj(repC));
                      case 'trivialRows'
                        repR = repC.group.trivialRep(repC.field, repC.dimension);
                      case 'trivialCols'
                        error('The row representation cannot be omitted for the trivialCols space');
                    end
                end
                if isempty(repC)
                    switch special
                      case 'antilinear'
                        repC = conj(repR);
                      case {'bilinear', 'symmetric'}
                        repC = dual(repR);
                      case 'commutant'
                        repC = repR;
                      case {'sesquilinear', 'hermitian'}
                        repC = dual(conj(repR));
                      case 'trivialRows'
                        error('The row representation cannot be omitted for the trivialCols space');
                      case 'trivialCols'
                        repC = repR.group.trivialRep(repR.field, repR.dimension);
                    end
                end
            end
            assert(repR.group == repC.group, 'Both representations must be defined over the same group');
            assert(repR.field == repC.field, 'Both representations must be defined over the same field');
        end

        function E = make(repR, repC, varargin)
        % Returns the space of equivariant linear maps between two representations
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`.Rep`): Representation on the target/row space
        %   repC (`.Rep`): Representation on the source/column space
        %
        % Keyword Args:
        %   special (charstring, optional): Special structure if applicable, see `.Equivariant`, default: ''
        %   type ('exact', 'double' or 'double/sparse', optional): Whether to obtain an exact equivariant space, default 'double' ('double' and 'double/sparse' are equivalent)
        %
        % Returns:
        %   `.Equivariant`: The equivariant vector space
            args = struct('special', '', 'type', 'double');
            args = replab.util.populateStruct(args, varargin);
            [repR, repC] = replab.Equivariant.validateArguments(repR, repC, args.special);
            if isa(repR.group, 'replab.FiniteGroup')
                E = replab.equi.Equivariant_forMonomialRep.make(repR, repC, args.special);
                if ~isempty(E)
                    return
                end
            end
            if isa(repR, 'replab.SubRep') && isa(repC, 'replab.SubRep')
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
