classdef NiceFiniteSubgroup < replab.NiceFiniteGroup
% A generic implementation of a subgroup of a nice finite group

    methods

        function self = NiceFiniteSubgroup(type, generators, order)
        % Constructs a subgroup of a nice finite group
        %
        % Args:
        %   type (`replab.NiceFiniteGroup`, optional): Type of this group, must not be a `NiceFiniteSubgroup`
        %   generators (cell(1,\*) of `.type` elements): Group generators
        %   order (vpi, optional): Order of the group
            self.type = type;
            self.representative = type.identity;
            self.identity = type.identity;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.cache('order', order, '==');
            end
            for i = 1:length(generators)
                assert(~type.isIdentity(generators{i}), 'Generator cannot be identity');
            end
            self.generators = generators;
        end

    end

    methods % Implementations

        % Domain

        function b = eqv(self, x, y)
            b = self.type.eqv(x, y);
        end

        % Monoid

        function z = compose(self, x, y)
            z = self.type.compose(x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = self.type.inverse(x);
        end

        % NiceFiniteGroup

        function p = niceImage(self, g)
            p = self.type.niceImage(g);
        end

    end

end
