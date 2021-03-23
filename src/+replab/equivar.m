classdef equivar < replab.Str
% Describes an equivariant map

    properties (SetAccess = protected)
        field % ('R', 'C'): Field
        repR % (`.Rep`): Row representation
        repC % (`.Rep`): Column representation
        special % ('', 'symmetric', 'hermitian'): Type of equivariant variable
        equivariant % (`.IrreducibleEquivariant`): Irreducible equivariant space
        blocks % (cell(\*,\*) of sdpvar(\*,\*,\*)): Matrix blocks
    end

    methods (Static)

        function [repR, repC] = validateReps(repR, repC, special)
        % Validate and complete constructor arguments
        %
        % Args:
        %   repR (`.Rep` or ``[]``): Representation acting on rows
        %   repC (`.Rep` or ``[]``): Representation acting on columns
        %   special ('', 'symmetric', 'hermitian'): Matrix type
        %
        % Returns
        % -------
        %   repR: `.Rep`
        %     Row representation
        %   repC: `.Rep`
        %     Column representation
            assert(ismember(special, {'' 'symmetric' 'hermitian'}), 'Special must be empty, ''symmetric'' or ''hermitian''.');
            assert(~isempty(repR) || ~isempty(repC), 'One of repR or repC must be provided');
            if isempty(repR)
                switch special
                  case 'symmetric'
                    assert(repC.overR, 'Symmetric equivariant matrices must be defined using real reps.');
                    repR = repC.';
                  case 'hermitian'
                    assert(repC.overC, 'Hermitian equivariant matrices must be defined using complex reps.');
                    repR = repC';
                  otherwise
                    error('If one rep argument is omitted, the matrix must be symmetric or hermitian');
                end
            end
            if isempty(repC)
                switch special
                  case 'symmetric'
                    assert(repR.overR, 'Symmetric equivariant matrices must be defined using real reps.');
                    repC = repR.';
                  case 'hermitian'
                    assert(repR.overC, 'Hermitian equivariant matrices must be defined using complex reps.');
                    repC = repR';
                  otherwise
                    error('If one rep argument is omitted, the matrix must be symmetric or hermitian');
                end
            end
            assert(repR.group == repC.group, 'Both representations must be defined over the same group');
            assert(repR.field == repC.field, 'Both representations must be defined over the same field');
        end

    end

    methods

        function self = equivar(repR, repC, special, equivariant, blocks)
        % Constructs an equivariant YALMIP matrix variable
        %
        % If a symmetric or hermitian matrix is requested, then one of the representation can be omitted and will
        % be deduced.
        %
        % Args:
        %   repR (`.Rep` or ``[]``): Representation acting on rows
        %   repC (`.Rep` or ``[]``): Representation acting on columns
        %   special ('', 'symmetric', 'hermitian'): Matrix type
        %   equivariant (`.IrreducibleEquivariant` or ``[]``): Equivariant space from ``repC.decomposition`` to ``repR.decomposition``, may be omitted
        %   blocks (cell(\*,\*) of sdpvar(\*,\*,\*)): Data present in each block
            if nargin < 3
                special = '';
            end
            if nargin < 4
                blocks = [];
            end
            if nargin < 5
                equivariant = [];
            end
            [repR, repC] = replab.equivar.validateReps(repR, repC, special);
            field = repR.field;
            if isempty(equivariant)
                switch special
                  case 'symmetric'
                    equivariant = repR.decomposition.sesquilinearInvariant;
                  case 'hermitian'
                    equivariant = repR.decomposition.sesquilinearInvariant;
                  case ''
                    equivariant = repR.decomposition.equivariantFrom(repC.decomposition);
                end
            end
            if isempty(blocks)
                B = equivariant.blocks;
                n1 = size(B, 1);
                n2 = size(B, 2);
                blocks = cell(n1, n2);
                for i = 1:size(B, 1)
                    for j = 1:size(B, 2)
                        if equivariant.nonZeroBlocks(i, j)
                            blocks{i,j} = equivariant.blocks{i,j}.makeSdpvar(special);
                        end
                    end
                end
            end
            self.field = field;
            self.repR = repR;
            self.repC = repC;
            self.special = special;
            self.equivariant = equivariant;
            self.blocks = blocks;
        end

        function C = sdp(self)
        % Returns the YALMIP constraint that this equivar is semidefinite positive
        %
        % This expands the SDP constraint in the block-diagonal basis
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
                    row{j} = self.equivariant.blocks{i,j}.reconstruct(self.blocks{i,j});
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
