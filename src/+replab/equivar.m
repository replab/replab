classdef equivar < replab.Obj
% Describes an equivariant YALMIP matrix variable
%
% The variable is defined over an equivariant space `.equivariant`. Internally, the variable will be parameterized using
% the irreducible decomposition of both ``equivariant.repR`` and ``equivariant.repC``, which is given by
% ``.equivariant.decomposition``, of type `.IrreducibleEquivariant`, which provides us a minimal parameterization
% of the different blocks.

    properties (SetAccess = protected)
        equivariant % (`.Equivariant`): Equivariant space corresponding to this matrix
        blocks % (cell(\*,\*) of double(\*,\*,\*) or sdpvar(\*,\*,\*)): Matrix blocks
    end

    methods (Static)

        function [lhs, rhs, E] = shapeArgs(lhs, rhs)
            if isa(lhs, 'replab.equivar') && isa(rhs, 'replab.equivar')
                assert(lhs.equivariant == rhs.equivariant);
                E = lhs.equivariant;
            elseif isa(lhs, 'replab.equivar')
                E = lhs.equivariant;
            elseif isa(rhs, 'replab.equivar')
                E = rhs.equivariant;
            end
            if isa(lhs, 'double')
                if isscalar(lhs) && lhs == 0
                    lhs = zeros(E.repR.dimension, E.repC.dimension);
                end
                lhs = replab.equivar(E, 'value', lhs);
            end
            if isa(rhs, 'double')
                if isscalar(rhs) && rhs == 0
                    rhs = zeros(E.repR.dimension, E.repC.dimension);
                end
                rhs = replab.equivar(E, 'value', rhs);
            end
            assert(isa(lhs, 'replab.equivar') && isa(rhs, 'replab.equivar'));
        end

    end

    methods

        function self = equivar(equivariant, varargin)
        % Constructs an equivariant YALMIP matrix variable
        %
        % There are two possibilities:
        %
        % * if one of the ``value`` or ``blocks`` arguments is given, this is initialized using the argument contents,
        % * if none of these keyword arguments is provided, the variable is free, and all its degrees of freedom will be
        %   initialized using YALMIP variables.
        %
        % Args:
        %   equivariant (`.Equivariant`): Equivariant space over which this variable is defined
        %
        % Keyword Args:
        %   value (double(\*,\*) or sdpvar(\*,\*)): Variable value (in the original space)
        %   blocks (cell(\*,\*) of double(\*,\*) or sdpvar(\*,\*)): Variable value in the minimal parameter space


        % TODO: cache "sdpvar(equivar)"
        % prefill that cache if value is provided here
        % allow an option to *not* project before factorization
            repR = equivariant.repR;
            repC = equivariant.repC;
            args = struct('value', [], 'blocks', []);
            args = replab.util.populateStruct(args, varargin);
            nProvided = ~isempty(args.value) + ~isempty(args.blocks);
            assert(nProvided <= 1, 'Provide at most one of the following keyword arguments: value, symmetryAdaptedValue, blocks');
            dec = equivariant.decomposition;
            if nProvided == 0
                blocks = dec.makeSdpvarBlocks;
            elseif ~isempty(args.blocks) % the value is known and in the internal format
                blocks = args.blocks;
            else
                n1 = dec.nRowBlocks;
                n2 = dec.nColBlocks;
                blocks = cell(n1, n2);
                if ~isempty(args.value)
                    self.cache('value', args.value, 'error');
                    value = args.value;
                    for i = 1:n1
                        for j = 1:n2
                            if isa(value, 'double')
                                [M, err] = dec.blocks{i,j}.projectAndFactorFromParent(value);
                            else
                                [M, err] = replab.evar.projectAndFactorFromParent_sdpvar(dec.blocks{i,j}, value);
                            end
                            blocks{i,j} = M;
                            replab.msg(2, 'equivar construction, block (%d,%d): maximum error %6.2f', i, j, max(err));
                        end
                    end
                end
            end
            self.equivariant = equivariant;
            self.blocks = blocks;
        end

        function rep = repR(self)
            rep = self.equivariant.repR;
        end

        function rep = repC(self)
            rep = self.equivariant.repC;
        end

        function s = double(self)
        % Returns the explicit double floating-point matrix corresponding to this equivar
        %
        % Raises:
        %   An error if the matrix is of type "sdpvar"
        %
        % Returns:
        %   double(\*,\*): Matrix
            s = self.sdpvar;
            assert(isa(s, 'double'));
        end

        function s = sdpvar(self)
        % Returns the explicit matrix corresponding to this equivar in the non-symmetry-adapted basis
        %
        % Returns:
        %   double(\*,\*) or sdpvar(\*,\*): Matrix
            s = self.cached('value', @() self.computeSdpvar);
        end

    end

    methods (Access = protected)

        function s = computeSdpvar(self)
            nR = size(self.blocks, 1);
            nC = size(self.blocks, 2);
            dec = self.equivariant.decomposition;
            values = cell(nR, 1);
            for i = 1:nR
                row = cell(1, nC);
                for j = 1:nC
                    row{j} = dec.blocks{i,j}.reconstruct(self.blocks{i,j});
                end
                values{i} = horzcat(row{:});
            end
            values = vertcat(values{:});
            s = self.equivariant.decomposition.repR.injection * values * self.equivariant.decomposition.repC.projection;
            switch self.equivariant.special
              case 'symmetric'
                s = (s + s.')/2;
              case 'hermitian'
                s = (s + s')/2;
            end
        end

    end

    methods % Arithmetic

        function res = plus(lhs, rhs)
            [lhs, rhs] = replab.equivar.shapeArgs(lhs, rhs);
            n1 = size(lhs.blocks, 1);
            n2 = size(lhs.blocks, 2);
            blocks = cell(n1, n2);
            for i = 1:n1
                for j = 1:n2
                    blocks{i,j} = lhs.blocks{i,j} + rhs.blocks{i,j};
                end
            end
            res = replab.equivar(lhs.equivariant, 'blocks', blocks);
        end

        function res = uminus(self)
            n1 = size(self.blocks, 1);
            n2 = size(self.blocks, 2);
            blocks = cell(n1, n2);
            for i = 1:n1
                for j = 1:n2
                    blocks{i,j} = -self.blocks{i,j};
                end
            end
            res = replab.equivar(self.equivariant, 'blocks', blocks);
        end

        function res = minus(lhs, rhs)
            res = lhs + (-rhs);
        end

        function res = mtimes(lhs, rhs)
            lhs_scalar = (isa(lhs, 'double') || isa(lhs, 'sdpvar')) && isscalar(lhs);
            rhs_scalar = (isa(rhs, 'double') || isa(rhs, 'sdpvar')) && isscalar(rhs);

            if isa(lhs, 'replab.equivar') && rhs_scalar
                n1 = size(lhs.blocks, 1);
                n2 = size(lhs.blocks, 2);
                blocks = cell(n1, n2);
                for i = 1:n1
                    for j = 1:n2
                        blocks{i,j} = lhs.blocks{i,j} * rhs;
                    end
                end
                res = replab.equivar(lhs.equivariant, 'blocks', blocks);
            elseif lhs_scalar && isa(rhs, 'replab.equivar')
                res = rhs * lhs;
            else
                error('Invalid syntax');
            end
        end

    end

    methods % Constraint construction

        function C = eq(lhs, rhs)
            [lhs, rhs, E] = replab.equivar.shapeArgs(lhs, rhs);
            dec = E.decomposition;
            n1 = dec.repR.nComponents;
            n2 = dec.repC.nComponents;
            C = {};
            for i = 1:n1
                for j = 1:n2
                    if dec.nonZeroBlock(i,j)
                        L = lhs.blocks{i,j};
                        R = rhs.blocks{i,j};
                        diff = L - R;
                        if isa(diff, 'double')
                            err = 0;
                            for i = 1:size(L, 3)
                                err = err + norm(diff(:,:,i), 'fro');
                            end
                            tol = replab.globals.doubleEigTol;
                            if err > tol
                                error('Equality constraint infeasible at construction error %6.2f  > tol %6.2f', err, tol);
                            end
                            replab.msg(1, 'Equality constraint involving constant terms is satisfied up to norm %6.2f', err);
                        else
                            for i = 1:size(L, 3)
                                C{1,end+1} = (L(:,:,i) == R(:,:,i));
                            end
                        end
                    end
                end
            end
            C = vertcat(C{:});
        end

        function C = issdp(self)
        % Returns the YALMIP constraint that this equivar is semidefinite positive
        %
        % This expands the SDP constraint in the block-diagonal basis
        %
        % Returns:
        %
            if self.equivariant.field == 'R'
                assert(ismember(self.equivariant.special, {'hermitian' 'symmetric'}));
            else
                assert(strcmp(self.equivariant.special, 'hermitian'));
            end
            dec = self.equivariant.decomposition;
            n1 = dec.repR.nComponents;
            n2 = dec.repC.nComponents;
            C = cell(1, n1);
            for i = 1:n1
                M = self.blocks{i,i};
                A = dec.blocks{i,i}.A;
                B = kron(M(:,:,1), A(:,:,1));
                for j = 2:size(M, 3)
                    B = B + kron(M(:,:,i), A(:,:,i));
                end
                B = (B + B')/2;
                C{i} = B >= 0;
            end
            C = vertcat(C{:});
        end

    end

end
