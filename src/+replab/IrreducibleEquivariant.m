classdef IrreducibleEquivariant < replab.SubEquivariant
% Equivariant space between two harmonized irreducible decompositions
%
% Matrices in this equivariant space are composed of blocks with the following form:
%
% $ X_ij = \sum_k M_ijk \otimes R_ijk \otimes A_ijk $ (or, most often, $X_ij$ is identically zero)
%
% where
%
% - $M_ijk$ represents the multiplicity space,
% - $R_ijk$ is a constant matrix representing the representation space,
% - $A_ijk$ is a constant matrix encoding the division algebra.

    properties (SetAccess = protected)
        blockRowSize % (integer(1,\*)): Row sizes of blocks
        blockColSize % (integer(1,\*)): Column sizes of blocks
        blocks % (cell(\*,\*) of `.IsotypicEquivariant` or ``[]``): Isotypic equivariant spaces
        nonZeroBlocks % (logical(\*,\*)): True when the corresponding block can be nonzero
    end


    methods (Access = protected)

        function self = IrreducibleEquivariant(parent, repR, repC, special, blocks)
            self@replab.SubEquivariant(parent, repR, repC, special);
            self.blocks = blocks;
            self.blockRowSize = cellfun(@(c) c.dimension, repR.components);
            self.blockColSize = cellfun(@(c) c.dimension, repC.components);
            self.nonZeroBlocks = cellfun(@(b) ~isempty(b), blocks);
        end

    end

    methods % Factorization

        function X = reconstruct(self, M, type)
        % Reconstructs an equivariant space element returned by a factorization method
        %
        % Args:
        %   M (cell(\*,\*) of cell(1,\*) of double(\*,\*) or cyclotomic(\*,\*)): Factorized element
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Equivariant matrix
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            rows = cell(1, self.repR.nComponents);
            for i = 1:self.repR.nComponents
                cols = cell(1, self.repC.nComponents);
                for j = 1:self.repC.nComponents
                    if self.nonZeroBlocks(i, j)
                        cols{j} = self.blocks{i, j}.reconstruct(M{i, j}, type);
                    else
                        if strcmp(type, 'exact')
                            cols{j} = replab.cyclotomic.zeros(self.blockRowSize(i), self.blockColSize(j));
                        else
                            cols{j} = zeros(self.blockRowSize(i), self.blockColSize(j));
                        end
                    end
                end
                rows{i} = horzcat(cols{:});
            end
            X = vertcat(rows{:});
        end

    end

    methods (Access = protected)

        function M = projectAndFactor_exact(self, X)
            M = cell(self.repR.nComponents, self.repC.nComponents);
            for i = 1:self.repR.nComponents
                for j = 1:self.repC.nComponents
                    if self.nonZeroBlocks(i, j)
                        rangeR = sum(self.blockRowSize(1:i-1)) + (1:self.blockRowSize(i));
                        rangeC = sum(self.blockColSize(1:j-1)) + (1:self.blockRowSize(j));
                        M{i, j} = self.blocks{i, j}.projectAndFactor(X(rangeR, rangeC), 'exact');
                    end
                end
            end
        end

        function M = projectAndFactor_double_sparse(self, X)
            M = cell(self.repR.nComponents, self.repC.nComponents);
            for i = 1:self.repR.nComponents
                for j = 1:self.repC.nComponents
                    if self.nonZeroBlocks(i, j)
                        rangeR = sum(self.blockRowSize(1:i-1)) + (1:self.blockRowSize(i));
                        rangeC = sum(self.blockColSize(1:j-1)) + (1:self.blockRowSize(j));
                        M{i, j} = self.blocks{i, j}.projectAndFactor(X(rangeR, rangeC), 'double/sparse');
                    end
                end
            end
        end

        function M = projectAndFactorFromParent_exact(self, parentX)
            M = cell(self.repR.nComponents, self.repC.nComponents);
            for i = 1:self.repR.nComponents
                for j = 1:self.repC.nComponents
                    if self.nonZeroBlocks(i, j)
                        M{i, j} = self.blocks{i, j}.projectAndFactorFromParent(X, 'exact');
                    end
                end
            end
        end

        function M = projectAndFactorFromParent_double_sparse(self, parentX)
            M = cell(self.repR.nComponents, self.repC.nComponents);
            for i = 1:self.repR.nComponents
                for j = 1:self.repC.nComponents
                    if self.nonZeroBlocks(i, j)
                        M{i, j} = self.blocks{i, j}.projectAndFactorFromParent(X, 'double/sparse');
                    end
                end
            end
        end

    end

    methods % Projection

        function [M, err] = projectAndFactor(self, X, type)
        % Projects the given matrix in the commutant algebra and factors it
        %
        % It returns the decomposition of the projection.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the isotypic component space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   M: cell(\*,\*) of (``[]``, or cell(1,\*) of (double(\*,\*) or `.cyclotomic`(\*,\*)))
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'exact'
                M = self.projectAndFactor_exact(X);
                err = 0;
              case 'double'
                M = self.projectAndFactor_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
              case 'double/sparse'
                M = self.projectAndFactor_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
            if nargout > 1
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X);
                err = sX*(eR*cC + cR*eC);
            end
        end

    end

    methods % Projection from parent space

        function [M, err] = projectAndFactorFromParent(self, X, type)
        % Projects the given matrix in the parent representation space into the commutant algebra and factors it
        %
        % It returns the decomposition of the projection.
        %
        % Args:
        %   X (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the parent representation space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   M: cell(\*,\*) of (``[]``, or cell(1,\*) of (double(\*,\*) or `.cyclotomic`(\*,\*)))
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
                case 'exact'
                  M = self.projectAndFactorFromParent_exact(X);
                  err = 0;
              case 'double'
                M = self.projectAndFactorFromParent_double_sparse(X);
                M = cellfun(@full, M, 'uniform', 0);
              case 'double/sparse'
                M = self.projectAndFactorFromParent_double_sparse(X);
              otherwise
                error('Unknown type %s', type);
            end
            if nargout > 1
                assert(replab.globals.yolo);
                eR = self.repR.errorBound;
                eC = self.repC.errorBound;
                cR = self.repR.conditionNumberEstimate; % condition number of repR
                cC = self.repC.conditionNumberEstimate; % condition number of repC
                sX = 0;
                R = self.R;
                A = self.A;
                X1 = kron(M{1}, kron(R{1}, A{1}));
                for i = 1:length(M)
                    sX = sX + replab.numerical.norm2UpperBound(M{i}) * replab.numerical.norm2UpperBound(R{i}) * replab.numerical.norm2UpperBound(A{i});
                end
                err = sX*(eR*cC + cR*eC);
            end
        end

    end


    methods (Static)

        function E = make_exact(parent, repR, repC, special)
            error('Not implemented'); % TODO
        end

        function E = make_double(parent, repR, repC, special, parentSample)
            if isempty(parentSample)
                parentSample = parent.sample;
            end
            blocks = cell(repR.nComponents, repC.nComponents);
            todo = true(repR.nComponents, repC.nComponents);
            for i = 1:repR.nComponents
                for j = 1:repC.nComponents
                    if todo(i, j)
                        E = replab.IsotypicEquivariant.make(repR.component(i), repC.component(j), '', 'parent', parent, 'parentSample', parentSample, 'type', 'double');
                        blocks{i, j} = E;
                        if ~isempty(E)
                            todo(i, :) = false;
                            todo(:, j) = false;
                        end
                    end
                end
            end
            E = replab.IrreducibleEquivariant(parent, repR, repC, special, blocks);
        end

        function E = make(repR, repC, special, varargin)
        % Returns the space of equivariant linear maps between two harmonized irreducible decompositions
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`+replab.Irreducible`): Irreducible decomposition acting on the target/row space
        %   repC (`+replab.Irreducible`): Irreducible decompositoin acting on the source/column space
        %   special (charstring): Special structure see help on `+replab.Equivariant.special`
        %
        % Keyword Args:
        %   type ('exact', 'double' or 'double/sparse'): Whether to obtain an exact equivariant space ('double' and 'double/sparse' are equivalent)
        %   parent (`.Equivariant`, optional): Equivariant space from ``repC.parent`` to ``repR.parent``
        %   parentSample (double(\*,\*), optional): A generic sample from the parent space
        %
        % Returns:
        %   `+replab.IrreducibleEquivariant`: The equivariant vector space
            assert(isa(repR, 'replab.Irreducible'));
            assert(isa(repC, 'replab.Irreducible'));
            assert(repR.isHarmonized);
            assert(repC.isHarmonized);
            assert(repR.field == repC.field);
            args = struct('parent', [], 'type', 'double', 'parentSample', []);
            args = replab.util.populateStruct(args, varargin);
            parent = args.parent;
            if isempty(parent)
                parent = repR.parent.equivariantFrom(repC.parent);
            end
            if strcmp(args.type, 'exact')
                assert(repR.isExact);
                assert(repC.isExact);
                assert(parent.isExact);
                E = replab.IrreducibleEquivariant.make_exact(parent, repR, repC, special);
            elseif strcmp(args.type, 'double') || strcmp(args.type, 'double/sparse')
                E = replab.IrreducibleEquivariant.make_double(parent, repR, repC, special, args.parentSample);
            else
                error('Wrong type value');
            end
        end

    end

end
