classdef equivar < replab.Str
% Describes an equivariant YALMIP matrix variable
%
% The variable is defined over an equivariant space `.equivariant`. Internally, the variable will be parameterized using
% the irreducible decomposition of both ``equivariant.repR`` and ``equivariant.repC``. Thus, we store the
% symmetry adapted equivariant space `.symmetryAdaptedEquivariant`, which is the equivariant space from
% ``equivariant.repC.decomposition`` to ``equivariant.repR.decomposition``.
%
% This symmetry adapted equivariant space is of type `.IrreducibleEquivariant`, which provides us a minimal parameterization
% of the different blocks.

    properties (SetAccess = protected)
        equivariant % (`.Equivariant`): Equivariant space corresponding to this matrix
        symmetryAdaptedEquivariant % (`.Equivariant`): Equivariant space in the symmetry adapted basis
        blocks % (cell(\*,\*) of double(\*,\*,\*) or sdpvar(\*,\*,\*)): Matrix blocks
    end

    methods

        function self = equivar(equivariant, varargin)
        % Constructs an equivariant YALMIP matrix variable
        %
        % There are two possibilities:
        %
        % * if one of the ``value``, ``symmetryAdaptedValue``, ``blocks`` keyword arguments is provided,
        %   the variable is initialized using the argument contents,
        % * if none of these keyword arguments is provided, the variable is free, and all its degrees of freedom will be
        %   initialized using YALMIP variables.
        %
        % If the `.symmetryAdaptedEquivariant` space has already been computed, it can be provided as well.
        %
        % Args:
        %   equivariant (`.Equivariant`): Equivariant space over which this variable is defined
        %
        % Keyword Args:
        %   value (double(\*,\*) or sdpvar(\*,\*)): Variable value (in the original space)
        %   symmetryAdaptedValue (double(\*,\*) or sdpvar(\*,\*)): Variable value (in the space of `.symmetryAdaptedEquivariant`)
        %   blocks (cell(\*,\*) of double(\*,\*) or sdpvar(\*,\*)): Variable value in the minimal parameter space
        %   symmetryAdaptedEquivariant (`.IrreducibleEquivariant`): Symmetry adapted equivariant space, see `.equivar`
            repR = equivariant.repR;
            repC = equivariant.repC;
            args = struct('value', [], 'symmetryAdaptedValue', [], 'blocks', [], 'symmetryAdaptedEquivariant', []);
            args = replab.util.populateStruct(args, varargin);
            symmetryAdaptedEquivariant = args.symmetryAdaptedEquivariant;
            if isempty(symmetryAdaptedEquivariant)
                switch equivariant.special
                  case 'antilinear'
                    symmetryAdaptedEquivariant = repR.decomposition.antilinearInvariant;
                  case 'bilinear'
                    symmetryAdaptedEquivariant = repC.decomposition.bilinearInvariant;
                  case 'commutant'
                    symmetryAdaptedEquivariant = repR.decomposition.commutant;
                  case 'hermitian'
                    symmetryAdaptedEquivariant = repC.decomposition.hermitianInvariant;
                  case 'sesquilinear'
                    symmetryAdaptedEquivariant = repC.decomposition.sesquilinearInvariant;
                  case 'symmetric'
                    symmetryAdaptedEquivariant = repC.decomposition.symmetricInvariant;
                  case 'trivialRows'
                    symmetryAdaptedEquivariant = repC.decomposition.trivialRowSpace;
                  case 'trivialCols'
                    symmetryAdaptedEquivariant = repR.decomposition.trivialColSpace;
                  case ''
                    symmetryAdaptedEquivariant = repR.decomposition.irreducibleEquivariantFrom(repC.decomposition);
                end
            else
                assert(strcmp(equivariant.special, symmetryAdaptedEquivariant.special));
            end
            nProvided = ~isempty(args.value) + ~isempty(args.symmetryAdaptedValue) + ~isempty(args.blocks);
            assert(nProvided <= 1, 'Provide at most one of the following keyword arguments: value, symmetryAdaptedValue, blocks');
            if nProvided == 0
                blocks = symmetryAdaptedEquivariant.makeSdpvarBlocks;
            elseif ~isempty(args.blocks) % the value is known and in the internal format
                blocks = args.blocks;
            else
                n1 = repR.decomposition.nComponents;
                n2 = repC.decomposition.nComponents;
                blocks = cell(n1, n2);
                if ~isempty(args.value)
                    value = args.value;
                    for i = 1:n1
                        for j = 1:n2
                            [M, err] = symmetryAdaptedEquivariant.blocks{i,j}.projectAndFactorFromParent(value);
                            blocks{i,j} = M;
                            err
                        end
                    end
                else % ~isempty(args.symmetryAdaptedValue)
                    symmetryAdaptedValue = args.symmetryAdaptedValue;
                    for i = 1:n1
                        for j = 1:n2
                            IE = symmetryAdaptedEquivariant.blocks{i,j};
                            shift1 = sum(cellfun(@(iso) iso.dimension, repR.decomposition.components(1:n1-1)));
                            shift2 = sum(cellfun(@(iso) iso.dimension, repC.decomposition.components(1:n2-1)));
                            range1 = shift1+(1:repR.decomposition.component(i).dimension);
                            range2 = shift2+(1:repC.decomposition.component(j).dimension);
                            B = symmetryAdaptedValue(range1, range2);
                            [M, err] = symmetryAdaptedEquivariant.projectAndFactor(symmetryAdaptedValue(range1, range2));
                            blocks{i,j} = M;
                            err
                        end
                    end
                end
            end

            self.equivariant = equivariant;
            self.symmetryAdaptedEquivariant = symmetryAdaptedEquivariant;
            self.blocks = blocks;
        end

        function rep = repR(self)
            rep = self.equivariant.repR;
        end

        function rep = repC(self)
            rep = self.equivariant.repC;
        end

        function C = sdp(self)
        % Returns the YALMIP constraint that this equivar is semidefinite positive
        %
        % This expands the SDP constraint in the block-diagonal basis
            if self.equivariant.field == 'R'
                assert(ismember(self.equivariant.special, {'symmetric' 'hermitian'}));
            else % self.equivariant.field == 'C'
                assert(strcmp(self.equivariant.special, 'hermitian'));
            end
            C = arrayfun(@(i) self.blocks{i,i} >= 0, 1:size(self.blocks, 1), 'uniform', 0);
            C = horzcat(C{:});
        end

        function s = sdpvar(self)
        % Returns the matrix corresponding to this equivar, in the original basis
            nR = size(self.blocks, 1);
            nC = size(self.blocks, 2);
            values = cell(nR, 1);
            for i = 1:nR
                row = cell(1, nC);
                for j = 1:nC
                    row{j} = self.symmetryAdaptedEquivariant.blocks{i,j}.reconstruct(self.blocks{i,j});
                end
                values{i} = horzcat(row{:});
            end
            values = vertcat(values{:});
            s = self.repR.decomposition.injection * values * self.repC.decomposition.projection;
        end

        function c = linearEqualityConstraint(self, F, Y)
        % Constructs a linear equality constraint taking in account the symmetries
        %
        % Returns a YALMIP linear equality constraint corresponding to the equation:
        % ``F(self.sdpvar) == Y.sdpvar``
        %
        % Args:
        %   F (function_handle): Linear map given as a generic function
        %   Y (`.equivar`): Right-hand side
        %
        % Returns:
        %   YALMIP constraint object: The constraint
            error('TODO');
        end

    end

end
