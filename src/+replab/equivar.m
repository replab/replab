classdef equivar < replab.Str
% Describes an equivariant map

    properties (SetAccess = protected)
        repR % (`.Rep`): Row representation
        repC % (`.Rep`): Column representation
        special % ('', 'symmetric', 'hermitian'): Type of equivariant variable
        blocks % (cell(1,\*) of sdpvar(\*,\*)): Matrix blocks
        equivariant % (`.IrreducibleEquivariant`): Irreducible equivariant space
    end

    methods

        function self = equivar(repR, repC, special, blocks, equivariant)
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
            if isempty(blocks)
                assert(~isempty(special), 'TODO: allow generic equivar');
            end
            self.repR = repR;
            self.repC = repC;
            self.special = special;
        end

    end

end
