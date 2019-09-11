classdef NiceFiniteSubgroup < replab.FiniteGroup
% A subgroup of a nice finite group
%
% Reuses much of the parent structure, including the group binary operation and the nice monomorphism.

    properties (SetAccess = protected)
        parent % replab.NiceFiniteGroup: Parent nice finite group
    end

    methods
        
        function self = NiceFiniteSubgroup(parent, generators, order)
        % Constructs a nice finite subgroup with the given properties
            self.parent = parent;
            self = self@replab.NiceFiniteGroup(parent.identity, generators, parent.niceMonomorphism);
        end
        
        % Domain
        
        function b = eqv(self, x, y)
            b = self.parent.eqv(x, y);
        end
        
        % Monoid
        
        function z = compose(self, x, y)
            z = self.parent.compose(x, y);
        end
a        
        % Group
        
        function xInv = inverse(self, x)
            xInv = self.parent.inverse(x);
        end
        
    end

end
