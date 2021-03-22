classdef equivar < replab.Str
% Describes an equivariant map

    properties (SetAccess = protected)
        field % ('R', 'C'): Field
        repR % (`.Rep`): Row representation
        repC % (`.Rep`): Column representation
        special % ('', 'symmetric', 'hermitian'): Type of equivariant variable
        equivariant % (`.IrreducibleEquivariant`): Irreducible equivariant space
        blocks % (cell(\*,\*) of sdpvar(\*,\*)): Matrix blocks
    end

    methods

        function self = equivar(repR, repC, special, equivariant, blocks)
            if nargin < 3
                special = '';
            end
            if nargin < 4
                blocks = [];
            end
            if nargin < 5
                equivariant = [];
            end
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
                    assert(repC.overC, 'Hermitian equivariant matrices must be defined using complex reps.');
                    repC = repR';
                  otherwise
                    error('If one rep argument is omitted, the matrix must be symmetric or hermitian');
                end
            end
            assert(repR.group == repC.group, 'Both representations must be defined over the same group');
            assert(repR.field == repC.field, 'Both representations must be defined over the same field');
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
                if isempty(special)
                    if field == 'R'
                        args = {'full'};
                    else
                        args = {'full' 'complex'};
                    end
                    for i = 1:size(B, 1)
                        for j = 1:size(B, 2)
                            b = B{i,j};
                            blocks{i,j} = sdpvar(b.repR.multiplicity, b.repC.multiplicity, args{:});
                        end
                    end
                else
                    if field == 'R'
                        args = {'symmetric'};
                    else
                        args = {'hermitian' 'complex'};
                    end
                    assert(n1 == n2);
                    assert(all(all(equivariant.nonZeroBlocks == logical(eye(n1)))));
                    for i = 1:n1
                        b = B{i,i};
                        m = b.repR.multiplicity;
                        blocks{i,i} = sdpvar(m, m, args{:});
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
            error('TODO');
        end

    end

end
