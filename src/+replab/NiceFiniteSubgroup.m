classdef NiceFiniteSubgroup < replab.NiceFiniteGroup
% A generic implementation of a subgroup of a nice finite group

    methods

        function self = NiceFiniteSubgroup(parent, generators, order)
        % Constructs a subgroup of a nice finite group
        %
        % Args:
        %   parent (`replab.PermutationGroup`, optional): Parent of this group, must not be a `NiceFiniteSubgroup`
        %   generators (cell(1,\*) of permutation): Group generators
        %   order (vpi, optional): Order of the group
            self.parent = parent;
            self.identity = parent.identity;
            % own stuff
            if nargin > 2 && ~isempty(order)
                self.cache('order', order, '==');
            end
            for i = 1:length(generators)
                assert(~parent.isIdentity(generators{i}), 'Generator cannot be identity');
            end
            self.generators = generators;
        end

        function res = hasSameTypeAs(self, rhs)
            if isa(rhs, 'replab.NiceFiniteSubgroup')
                res = self.parent.hasSameTypeAs(rhs.parent);
            else
                res = self.parent.hasSameTypeAs(rhs);
            end
        end

        function p = niceMonomorphismImage(self, g)
            p = self.parent.niceMonomorphismImage(g);
        end

        %% Domain methods

        function b = eqv(self, x, y)
            b = self.parent.eqv(x, y);
        end

        %% Monoid methods

        function z = compose(self, x, y)
            z = self.parent.compose(x, y);
        end

        %% Group methods

        function xInv = inverse(self, x)
            xInv = self.parent.inverse(x);
        end

    end

end
