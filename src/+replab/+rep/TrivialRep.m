classdef TrivialRep < replab.Rep
% Describes d copies of the real or complex trivial representation of a group

    properties
        identity % Stored copy of the identity matrix
    end

    methods

        function self = TrivialRep(group, field, dimension)
            assert(isa(group, 'replab.CompactGroup'));
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            self.isUnitary = true;
            self.identity = eye(self.dimension);
        end

        %% Rep methods

        function rho = image(self, g)
            rho = self.identity;
        end

        function rho = inverseImage(self, g)
            rho = self.identity;
        end

        function rho = sample(self)
            rho = self.identity;
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

        function rep = conj(self)
            rep = self;
        end

        function rep = dual(self)
            rep = self;
        end

    end

end
