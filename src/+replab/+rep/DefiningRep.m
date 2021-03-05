classdef DefiningRep < replab.Rep
% Defining representation of a `+replab.ClassicalCompactGroup`

    methods

        function self = DefiningRep(group, field)
            switch [group.algebra '/' field]
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
              otherwise
                error('Invalid');
            end
            self@replab.Rep(group, field, d, 'isUnitary', true, 'trivialDimension', 0, 'isIrreducible', true, 'divisionAlgebraName', divisionAlgebraName);
        end

    end

    methods (Access = protected) % Implementations

        function e = computeErrorBound(self)
            e = inf;
        end

        function rho = image_double_sparse(self, g)
            switch [self.group.algebra '/' self.field]
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

        function b = hasTorusImage(self)
            b = true;
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            n = self.group.n;
            W = [1 1i; 1 -1i]/sqrt(2);
            switch [self.group.algebra '/' self.field]
              case {'R/R', 'R/C'}
                if mod(n, 2) == 0
                    torusMap = kron(eye(n/2), [1;-1]);
                    torusInjection = kron(speye(n/2), W');
                    torusProjection = kron(speye(n/2), W);
                else
                    torusMap = [kron(eye((n-1)/2), [1;-1]); zeros(1, (n-1)/2)];
                    torusInjection = blkdiag(kron(speye((n-1)/2), W'), 1);
                    torusProjection = blkdiag(kron(speye((n-1)/2), W), 1);
                end
              case {'C/C', 'C/R'}
                if self.group.isSpecial
                    torusMap = [eye(n-1); -ones(1,n-1)];
                else
                    torusMap = eye(n);
                end
                d = size(torusMap, 1);
                if self.field == 'C'
                    torusInjection = speye(d);
                    torusProjection = speye(d);
                else
                    torusMap = kron(torusMap, [1;-1]);
                    torusInjection = kron(speye(d), W');
                    torusProjection = kron(speye(d), W);
                end
              case 'H/C'
                torusMap = kron(eye(n), [1;-1]);
                d = size(torusMap, 1);
                torusInjection = speye(d);
                torusProjection = speye(d);
            end
        end

    end

end
