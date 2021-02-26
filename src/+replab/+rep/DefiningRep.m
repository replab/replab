classdef DefiningRep < replab.Rep
% Defining representation of a unitary or orthogonal group

    properties (SetAccess = protected)
        fromDivisionRing % ('R', 'C', 'H'): Division ring over which the group element matrices are defined
    end

    methods

        function self = DefiningRep(group, field)
            switch class(group)
              case 'replab.OrthogonalGroup'
                fromDivisionRing = 'R';
              case 'replab.UnitaryGroup'
                fromDivisionRing = 'C';
              case 'replab.CompactSymplecticGroup'
                fromDivisionRing = 'H';
              otherwise
                error('Unknown group type %s', class(group));
            end
            switch [fromDivisionRing '/' field]
              case {'R/R', 'C/C', 'R/C'}
                d = group.n;
                divisionAlgebraName = '';
              case 'H/C'
                d = group.n*2;
                divisionAlgebraName = '';
              case 'H/R'
                d = group.n*4;
                divisionAlgebraName = 'quaternion.rep';
              case 'C/R'
                d = group.n*2;
                divisionAlgebraName = 'complex';
            end
            self@replab.Rep(group, field, d, 'isUnitary', true, 'trivialDimension', 0, 'isIrreducible', true);
            self.fromDivisionRing = fromDivisionRing;
        end

    end

    methods (Access = protected) % Implementations

        function e = computeErrorBound(self)
            e = inf;
        end

        function rho = image_double_sparse(self, g)
            switch [self.fromDivisionRing '/' self.field]
              case {'R/R', 'C/C', 'R/C'}
                rho = g;
              case 'H/C'
                rho = kron(g.X, [1 0; 0 0]) + ...
                      kron(conj(g.X), [0 0; 0 1]) + ...
                      kron(g.Y, [0 1; 0 0]) + ...
                      kron(-conj(g.Y), [0 0; 1 0]);
              case 'H/R'
                rho = kron(g.part1, eye(4)) + ...
                      kron(g.parti, [0 -1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 1 0]) + ...
                      kron(g.partj, [0 0 -1 0; 0 0 0 1; 1 0 0 0; 0 -1 0 0]) + ...
                      kron(g.partk, [0 0 0 -1; 0 0 -1 0; 0 1 0 0; 1 0 0 0]);
              case 'C/R'
                rho = kron(real(g), eye(2)) + ...
                      kron(imag(g), [0 -1; 1 0]);
            end
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false; % Note: redundant, as parent is already false
        end

        function p = invariantBlocks(self)
            p = replab.Partition.fromBlocks({1:self.dimension});
        end

        function b = hasMaximalTorusExponents(self)
            b = false;
            switch [self.fromDivisionRing '/' self.field]
              case 'C/C'
                b = true;
              case 'H/C'
                b = true;
            end
        end

        function [powers, blockIndex] = maximalTorusExponents(self)
            assert(self.hasMaximalTorusExponents);
            switch [self.fromDivisionRing '/' self.field]
              case 'C/C'
                powers = eye(self.dimension);
                blockIndex = 1:self.dimension;
              case 'H/C'
                powers = kron(eye(self.dimension/2), [1;-1]);
                blockIndex = 1:self.dimension;
            end
        end

    end

end
