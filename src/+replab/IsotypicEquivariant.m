classdef IsotypicEquivariant < replab.SubEquivariant
% Equivariant space between two harmonized isotypic components containing equivalent representations
%
% We consider two cases.
%
% If the two isotypic components contain equivalent irreducible representations, then `.isZero` is false and
% the matrices in this equivariant space have the following form:
%
% $ X = \sum_i M(:,:,i) \otimes R(:,:,i) \otimes A(:,:,i) $
%
% where
%
% - $M$ represents the multiplicity space,
% - $R$ is a constant matrix representing the representation space,
% - $A$ is a constant matrix encoding the division algebra.
%
% The division algebra matrix $A$ is given by:
%
% - ``A = 1`` if `.divisionAlgebraName` is ``''``
% - ``A(:,:,1) = [1 0; 0 1]``, ``A(:,:,2) = [0 -1; 1 0]`` if `.divisionAlgebraName` is ``'C->R'``
% - ``A(:,:,1) = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]``, ``A(:,:,2) = [0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0]``,
%   ``A(:,:,3) = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]``, ``A(:,:,4) = [0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0]``
%   if `.divisionAlgebraName` is ``'H->R:equivariant'``.
%
% If the two isotypic components contain inequivalent irreps, then `.isZero` is true. `.divisionAlgebraName` is ``'0'``,
% and ``R`` has size ``[m n 0]``.
%
% Example:
%   >>> S3 = replab.S(3);
%   >>> natDec = S3.naturalRep.decomposition;
%   >>> std1 = natDec.component(2);
%   >>> nat2Dec = S3.naturalRep.tensorPower(2).decomposition;
%   >>> std2 = nat2Dec.component(3);
%   >>> E = std1.isotypicEquivariantFrom(std2);
%   >>> X = E.sample;
%   >>> M = E.projectAndFactor(X);
%   >>> length(M)
%       1
%   >>> R = E.R('double');
%   >>> A = E.A('double');
%   >>> tol = 1e-14;
%   >>> norm(kron(M{1}, kron(R{1}, A{1})) - X, 'fro') <= tol
%       1
%   >>> norm(E.reconstruct(M) - X, 'fro') <= tol
%       1

    properties (SetAccess = protected)
        R_internal % (double(\*,\*,\*) or `.cyclotomic`(\*,\*,\*)): Representation space basis
        divisionAlgebraName % ('', 'C->R', 'H->R:equivariant'): Division algebra
    end

    methods

        function self = IsotypicEquivariant(parent, repR, repC, special, R_internal, divisionAlgebraName)
            assert(isa(repR, 'replab.Isotypic'));
            assert(isa(repC, 'replab.Isotypic'));
            switch divisionAlgebraName
              case '0'
                d = 0; % division algebra dimension
                ms = 1;
              case ''
                d = 1;
                ms = 1;
              case 'C->R'
                assert(repR.overR);
                d = 2;
                ms = 2;
              case 'H->R:equivariant'
                assert(repR.overR);
                d = 4;
                ms = 4;
            end
            r1 = repR.irrepDimension / ms;
            r2 = repC.irrepDimension / ms;
            assert(size(R_internal, 1) == r1);
            assert(size(R_internal, 2) == r2);
            assert(size(R_internal, 3) == d);
            assert(ismember(divisionAlgebraName, {'0' '' 'C->R' 'H->R:equivariant'}));
            self@replab.SubEquivariant(parent, repR, repC, special);
            self.R_internal = R_internal;
            self.divisionAlgebraName = divisionAlgebraName;
        end

    end

    methods % Implementations

        % Domain

        function l = laws(self)
            l = replab.laws.IsotypicEquivariantLaws(self);
        end

    end

    methods % Properties

        function d = divisionAlgebraBlockSize(self)
        % Returns the dimension of the matrix block encoding the division algebra
        %
        % As a special case, if the block is zero, the returned value is one.
        %
        % Returns:
        %   integer: Dimension of the matrix block
            switch self.divisionAlgebraName
              case '0'
                d = 1;
              case ''
                d = 1;
              case 'C->R'
                d = 2;
              case 'H->R:equivariant'
                d = 4;
            end
        end

        function d = divisionAlgebraDimension(self)
        % Returns the dimension of the division algebra encoded in this block
        %
        % As a special case, if the block is zero, the returned value is zero.
        %
        % Returns:
        %   integer: Dimension of the division algebra
            switch self.divisionAlgebraName
              case '0'
                d = 0;
              case ''
                d = 1;
              case 'C->R'
                d = 2;
              case 'H->R:equivariant'
                d = 4;
            end
        end

        function b = isZero(self)
        % Returns whether this equivariant space contains only the zero matrix
        %
        % This happens when the isotypic components `.repR` and `.repC` correspond to inequivalent irreducible representations
        %
        % Example:
        %   >>> S3 = replab.S(3);
        %   >>> rep = S3.naturalRep;
        %   >>> triv = rep.decomposition.component(1);
        %   >>> std = rep.decomposition.component(2);
        %   >>> E = triv.isotypicEquivariantFrom(std);
        %   >>> E.isZero
        %       1
        %
        % Returns:
        %   logical: True if the equivariant space is trivial
            b = strcmp(self.divisionAlgebraName, '0');
        end

    end

    methods % Factorization

        function res = R(self, type)
        % Returns the representation space basis
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*,\*) or `.cyclotomic`(\*,\*,\*): The representation space basis
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            res = replab.numerical.convert(self.R_internal, type);
        end

        function res = A(self, type)
        % Returns the division algebra basis
        %
        % Args:
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %  double(\*,\*) or `.cyclotomic`(\*,\*): The division algebra basis
            if nargin < 2 || isempty(type)
                type = 'double';
            end
            switch self.divisionAlgebraName
              case '0'
                A = zeros(1, 1, 0);
              case ''
                A = ones(1, 1, 1);
              case 'C->R'
                A = zeros(2, 2, 2);
                A(:,:,1) = [1 0; 0 1];
                A(:,:,2) = [0 -1; 1 0];
              case 'H->R:equivariant'
                A(:,:,1) = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
                A(:,:,2) = [0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0];
                A(:,:,3) = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0];
                A(:,:,4) = [0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0];
            end
            res = replab.numerical.convert(A, type);
        end

        function X = reconstruct(self, M, type)
        % Reconstructs an equivariant space element returned by a factorization method
        %
        % Args:
        %   M (double(\*,\*,\*) or cyclotomic(\*,\*,\*)): Factorized element
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns:
        %   double(\*,\*) or cyclotomic(\*,\*): Equivariant matrix
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if self.isZero
                X = zeros(self.repR.dimension, self.repC.dimension);
                if strcmp(type, 'exact')
                    X = replab.cyclotomic(X);
                end
            else
                R = self.R(type);
                A = self.A(type);
                X = kron(M(:,:,1), kron(R(:,:,1), A(:,:,1)));
                for i = 2:size(M, 3)
                    X = X + kron(M(:,:,1), kron(R(:,:,1), A(:,:,1)));
                end
            end
        end

    end

    methods % YALMIP helpers

        function M = makeSdpvar(self)
        % Returns a parameterization of the multiplicity space using YALMIP variables
        %
        % Returns:
        %   sdpvar(\*,\*,\*): Parameterization
            d1 = self.repR.multiplicity;
            d2 = self.repC.multiplicity;
            da = self.divisionAlgebraDimension;
            if self.isZero
                M = zeros(d1, d2, 0);
                return
            end
            if self.field == 'R'
                switch self.special
                  case {'hermitian', 'symmetric'}
                    M = sdpvar(d1, d2, 1, 'symmetric');
                    switch self.divisionAlgebraName
                      case ''
                      case 'C->R'
                        M = reshape(M, [d1*d2, 1]);
                        M = horzcat(M, zeros(d1*d2, 1));
                        M = reshape(M, [d1, d2, 2]);
                      case 'H->R:equivariant'
                        M = reshape(M, [d1*d2, 1]);
                        M = horzcat(M, zeros(d1*d2, 3));
                        M = reshape(M, [d1, d2, 4]);
                    end
                  otherwise
                    M = sdpvar(d1, d2, self.divisionAlgebraDimension, 'full');
                end
            else % self.field == 'C'
                switch self.special
                  case 'hermitian'
                    M = sdpvar(d1, d2, 'hermitian', 'complex');
                    M = reshape(M, [d1 d2 1]);
                  case 'symmetric'
                    warning('Is that case correct? please be extra careful');
                    M = sdpvar(d1, d2, 'symmetric', 'complex');
                    M = reshape(M, [d1 d2 1]);
                  otherwise
                    M = sdpvar(d1, d2, 1, 'full', 'complex');
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
        %   M: double(\*,\*,\*) or `.cyclotomic`(\*,\*,\*)
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if self.isZero
                M = zeros(self.repR.multiplicity, self.repC.multiplicity, 0);
                if strcmp(type, 'exact')
                    M = replab.cyclotomic(M);
                end
            else
                switch type
                  case 'exact'
                    parentX = self.repR.injection('exact') * X * self.repC.projection('exact');
                    M = self.projectAndFactorFromParent(parentX, 'exact');
                    err = 0;
                  case {'double', 'double/sparse'}
                    parentX = full(self.repR.injection('double/sparse') * X * self.repC.projection('double/sparse'));
                    M = self.projectAndFactorFromParent(parentX, 'double');
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

    end

    methods % Projection from parent space

        function [M, err] = projectAndFactorFromParent(self, parentX, type)
        % Projects the given matrix in the parent representation space into the commutant algebra and factors it
        %
        % It returns the decomposition of the projection.
        %
        % Args:
        %   parentX (double(\*,\*) or `.cyclotomic`(\*,\*), may be sparse): Matrix in the parent representation space to project
        %   type ('double', 'double/sparse' or 'exact', optional): Type of the returned value, default: 'double'
        %
        % Returns
        % -------
        %   M: double(\*,\*,\*) or `.cyclotomic`(\*,\*,\*)
        %     The part containing the degrees of freedom of the commutant algebra
        %   err: double
        %     Estimation of the numerical error, expressed as the distance of the returned projection to the invariant subspace in Frobenius norm
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            if self.isZero
                M = zeros(self.repR.multiplicity, self.repC.multiplicity, 0);
                if strcmp(type, 'exact')
                    M = replab.cyclotomic(M);
                end
                err = 0;
            else
                X = self.projectFromParent(parentX, type);
                R = self.R(type);
                M = zeros(self.repR.multiplicity, self.repC.multiplicity, self.divisionAlgebraDimension);
                if strcmp(type, 'exact')
                    M = replab.cyclotomic(M);
                end
                if self.field == 'C' || self.repR.irrep(1).frobeniusSchurIndicator == 1
                    M(:,:,1) = replab.IsotypicEquivariant.kronFirstFactor(X, R(:,:,1));
                elseif repR.irrep(1).frobeniusSchurIndicator == -2 % quaternion-type
                    [X1, X2, X3, X4] = replab.domain.QuaternionTypeMatrices.fromMatrix(X, 'commutant');
                    M(:,:,1) = replab.IsotypicEquivariant.kronFirstFactor(X1, R(:,:,1));
                    M(:,:,2) = replab.IsotypicEquivariant.kronFirstFactor(X2, R(:,:,2));
                    M(:,:,3) = replab.IsotypicEquivariant.kronFirstFactor(X3, R(:,:,3));
                    M(:,:,4) = replab.IsotypicEquivariant.kronFirstFactor(X4, R(:,:,4));
                elseif repR.irrep(1).frobeniusSchurIndicator == 0 % complex-type
                    [X1, X2] = replab.domain.ComplexTypeMatrices.fromMatrix(X);
                    M(:,:,1) = replab.IsotypicEquivariant.kronFirstFactor(X1, R(:,:,1));
                    M(:,:,2) = replab.IsotypicEquivariant.kronFirstFactor(X2, R(:,:,2));
                end
                if nargout > 1
                    if strcmp(type, 'exact')
                        err = 0;
                    else
                        assert(replab.globals.yolo);
                        eR = self.repR.errorBound;
                        eC = self.repC.errorBound;
                        cR = self.repR.conditionNumberEstimate; % condition number of repR
                        cC = self.repC.conditionNumberEstimate; % condition number of repC
                        sX = 0;
                        R = self.R;
                        A = self.A;
                        for i = 1:self.divisionAlgebraDimension
                            sX = sX + replab.numerical.norm2UpperBound(M(:,:,i)) * replab.numerical.norm2UpperBound(R(:,:,i)) * replab.numerical.norm2UpperBound(A(:,:,i));
                        end
                        err = sX*(eR*cC + cR*eC);
                    end
                end
            end
        end

    end

    methods (Static, Access = protected)

        function E = make_exact(parent, repR, repC, special)
            error('Not implemented'); % TODO
        end

        function E = make_double(parent, repR, repC, special, parentSample)
            if isempty(parentSample)
                sub = replab.SubEquivariant(parent, repR, repC, special);
                [X, err] = sub.sample; % TODO: proper error estimation here
            else
                X = repR.projection('double/sparse') * parentSample * repC.injection('double/sparse');
                assert(replab.globals.yolo);
                eR = repR.errorBound;
                eC = repC.errorBound;
                cR = repR.conditionNumberEstimate; % condition number of repR
                cC = repC.conditionNumberEstimate; % condition number of repC
                sX = replab.numerical.norm2UpperBound(X);
                err = sX*(eR*cC + cR*eC);
            end
            err = max(err, 10*sqrt(repR.dimension*repC.dimension)*eps);
            assert(err < 1e-10, 'Error to big');
            if norm(X, 'fro') < err
                E = [];
                return
            end
            if repR.overC || repR.irrep(1).frobeniusSchurIndicator == 1
                ird = repR.irrepDimension;
                R = replab.IsotypicEquivariant.kronSecondFactor(X, ird, ird);
                R = R/(trace(R)/ird);
                divisionAlgebraName = '';
                switch special
                  case 'symmetric'
                    R = (R+R.')/2;
                  case 'hermitian'
                    R = (R+R')/2;
                end
                R = reshape(R, [size(R, 1) size(R, 2) 1]);
            elseif repR.irrep(1).frobeniusSchurIndicator == -2 % quaternion-type

                % TODO: symmetric/hermitian support
                d = repR.irrepDimension/4;
                [X1, X2, X3, X4] = replab.domain.QuaternionTypeMatrices.fromMatrix(X, 'commutant');
                R1 = replab.IsotypicEquivariant.kronSecondFactor(X1, d, d);
                R2 = replab.IsotypicEquivariant.kronSecondFactor(X2, d, d);
                R3 = replab.IsotypicEquivariant.kronSecondFactor(X3, d, d);
                R4 = replab.IsotypicEquivariant.kronSecondFactor(X4, d, d);
                R1 = R1/norm(R1,'fro')*sqrt(ird);
                R2 = R2/norm(R2,'fro')*sqrt(ird);
                R3 = R3/norm(R3,'fro')*sqrt(ird);
                R4 = R4/norm(R4,'fro')*sqrt(ird);
                R = zeros(size(R1, 1), size(R1, 2), 4);
                R(:,:,1) = R1;
                R(:,:,2) = R2;
                R(:,:,3) = R3;
                R(:,:,4) = R4;
                divisionAlgebraName = 'H->R:equivariant';
            elseif repR.irrep(1).frobeniusSchurIndicator == 0 % complex-type
                % TODO: symmetric/hermitian support
                d = repR.irrepDimension/2;
                [X1, X2] = replab.domain.ComplexTypeMatrices.fromMatrix(X);
                R1 = replab.IsotypicEquivariant.kronSecondFactor(X1, d, d);
                R2 = replab.IsotypicEquivariant.kronSecondFactor(X2, d, d);
                R1 = R1/norm(R1,'fro')*sqrt(d);
                R2 = R2/norm(R2,'fro')*sqrt(d);
                R = zeros(size(R1, 1), size(R1, 2), 2);
                R(:,:,1) = R1;
                R(:,:,2) = R2;
                divisionAlgebraName = 'C->R';
            else
                error('Unknown type');
            end
            E = replab.IsotypicEquivariant(parent, repR, repC, special, R, divisionAlgebraName);
        end

    end

    methods (Static)

        function X1 = kronFirstFactor(X, X2)
        % Factors a Kronecker product, second step
        %
        % Assuming that X is approximately ``X = kron(X1, X2)``, with ``X2`` known,
        % this estimates the first factor ``X1`` and returns it.
        %
        % Args:
        %   X (double(\*,\*)): Matrix to factor
        %   X2 (double(\*,\*)): Second factor
        %
        % Returns:
        %   double(\*,\*): An estimation of the first factor in the Kronecker product
            F = sum(sum(conj(X2).*X2));
            nRows2 = size(X2, 1);
            nCols2 = size(X2, 2);
            nRows1 = size(X, 1) / nRows2;
            nCols1 = size(X, 2) / nCols2;
            X1 = zeros(nRows1, nCols1);
            for r = 1:nRows1
                Ir = (r-1)*nRows2+(1:nRows2);
                for c = 1:nCols1
                    Ic = (c-1)*nCols2+(1:nCols2);
                    M = X(Ir, Ic);
                    X1(r, c) = sum(sum(M.*conj(X2)))/F;
                end
            end
        end

        function X2 = kronSecondFactor(X, nRows2, nCols2)
        % Factors a Kronecker product, first step
        %
        % Assuming that X is approximately ``X = kron(X1, X2)``, with ``X2`` a matrix of size ``nRows2`` x ``nCols2``,
        % this estimates the factor ``X2`` and returns it.
        %
        % Args:
        %   X (double(\*,\*)): Matrix to factor
        %   nRows2 (integer): Number of rows of the second factor
        %   nCols2 (integer): Number of columns of the second factor
        %
        % Returns:
        %   double(\*,\*): An estimation of the second factor in the Kronecker product
            nRows1 = size(X, 1) / nRows2;
            nCols1 = size(X, 2) / nCols2;
            assert(round(nRows1) == nRows1 && round(nCols1) == nCols1);
            X2 = zeros(nRows2, nCols2);
            for r = 1:nRows1
                Ir = (r-1)*nRows2+(1:nRows2);
                for c = 1:nCols1
                    Ic = (c-1)*nCols2+(1:nCols2);
                    M = X(Ir, Ic);
                    X2 = X2 + M/sign(M(1,1));
                end
            end
        end

        function E = make(repR, repC, varargin)
        % Returns the space of equivariant linear maps between two isotypic components
        %
        % The equivariant vector space contains the matrices X such that
        %
        % ``repC.image(g) * X = X * repR.image(g)``
        %
        % Args:
        %   repR (`+replab.Isotypic`): Isotypic component on the target/row space
        %   repC (`+replab.Isotypic`): Isotypic component on the source/column space
        %
        % Keyword Args:
        %   special (charstring): Special structure, see help on `.Equivariant`
        %   type ('exact', 'double' or 'double/sparse'): Whether to obtain an exact equivariant space ('double' and 'double/sparse' are equivalent)
        %   parent (`.Equivariant`, optional): Equivariant space from ``repC.parent`` to ``repR.parent``
        %   parentSample (double(\*,\*), optional): A generic sample from the parent space
        %
        % Returns:
        %   `+replab.IsotypicEquivariant`: The equivariant vector space
            assert(isa(repR, 'replab.Isotypic'));
            assert(isa(repC, 'replab.Isotypic'));
            assert(repR.field == repC.field);
            args = struct('special', '', 'parent', [], 'type', 'double', 'parentSample', []);
            args = replab.util.populateStruct(args, varargin);
            parent = args.parent;
            if isempty(parent)
                switch args.special
                  case 'antilinear'
                    parent = repR.parent.antilinearInvariant(args.type);
                  case 'bilinear'
                    parent = repC.parent.bilinearInvariant(args.type);
                  case 'commutant'
                    parent = repR.parent.commutant(args.type);
                  case 'hermitian'
                    parent = repC.parent.hermitianInvariant(args.type);
                  case 'sesquilinear'
                    parent = repC.parent.sesquilinearInvariant(args.type);
                  case 'symmetric'
                    parent = repC.parent.symmetricInvariant(args.type);
                  case {'', 'trivialRows', 'trivialCols'}
                    parent = repR.parent.equivariantFrom(repC.parent, 'type', args.type);
                  otherwise
                    error('Invalid special structure %s', args.special);
                end
            end
            if repR.dimension == 0 || repC.dimension == 0
                E = [];
            elseif repR.irrepDimension ~= repC.irrepDimension
                E = [];
            elseif repR.irrep(1).frobeniusSchurIndicator ~= repC.irrep(1).frobeniusSchurIndicator
                E = [];
            else
                if strcmp(args.type, 'exact')
                    assert(repR.isExact);
                    assert(repC.isExact);
                    assert(parent.isExact);
                    E = replab.IsotypicEquivariant.make_exact(parent, repR, repC, args.special);
                elseif strcmp(args.type, 'double') || strcmp(args.type, 'double/sparse')
                    E = replab.IsotypicEquivariant.make_double(parent, repR, repC, args.special, args.parentSample);
                else
                    error('Wrong type value');
                end
            end
            if isempty(E)
                E = replab.IsotypicEquivariant(parent, repR, repC, args.special, zeros(repR.irrepDimension, repC.irrepDimension, 0), '0');
            end
        end

    end

end
