classdef TrivialRep < replab.Rep
% Describes d copies of the real or complex trivial representation of a group

    methods

        function self = TrivialRep(group, field, dimension)
            assert(isa(group, 'replab.CompactGroup'));
            self@replab.Rep(group, field, dimension, 'isUnitary', true, 'isIrreducible', dimension == 1, ...
                            'trivialDimension', dimension, 'frobeniusSchurIndicator', dimension, ...
                            'isDivisionAlgebraCanonical', strcmp(field, 'R') & dimension == 1);
        end

    end

    methods % Implementations

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        function b = canComputeType(self, type)
            b = true;
        end

        function rho = image(self, g, type)
            if nargin < 3 || isempty(type)
                type = 'double';
            end
            switch type
              case 'double'
                rho = eye(self.dimension);
              case 'intval'
                rho = intval(eye(self.dimension));
              case 'cyclotomic'
                rho = replab.cyclotomic.eye(self.dimension);
            end
        end

        function rho = inverseImage(self, g, varargin)
            rho = self.image(g, varargin{:});
        end

        function M = matrixRowAction(self, g, M)
        % do nothing to M
        end

        function M = matrixColAction(self, g, M)
        % do nothing to M
        end

        function complexRep = complexification(self)
            assert(self.overR, 'Representation should be real to start with');
            complexRep = replab.rep.TrivialRep(self.group, 'C', self.dimension);
        end

    end

end
