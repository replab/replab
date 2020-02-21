classdef TrivialRep < replab.Rep
% Describes d copies of the real or complex trivial representation of a group

    methods

        function self = TrivialRep(group, field, dimension)
            assert(isa(group, 'replab.CompactGroup'));
            % replab.Rep, immutable
            self.group = group;
            self.field = field;
            self.dimension = dimension;
            % replab.Rep, mutable
            self.isIrreducible = (self.dimension == 1);
            self.isUnitary = true;
            self.trivialDimension = dimension;
        end

        %% Str

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        %% Rep methods

        function rho = image_internal(self, g)
            rho = speye(self.dimension);
        end

        function rho = inverseImage_internal(self, g)
            rho = speye(self.dimension);
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
