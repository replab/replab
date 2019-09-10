classdef FiniteSubgroup < replab.FiniteGroup
% Finite groups that can construct subgroups and test for membership
        
    properties (SetAccess = protected)
        parent % replab.FiniteSubgroup: Parent group that is its own parent
    end

    methods
    
        function b = contains(self, g)
        % Tests whether this group contains the given parent group element
        %
        % Abstract in `replab.FiniteSubgroup`
        %
        % Args:
        %   g (element of `self.parent`): Element to test membership of
        %
        % Returns:
        %   logical: True if this group contains `g` and false otherwise
            error('Not implemented');
        end
        
        function sub = subgroup(self, generators, order)
        % Constructs a subgroup of the current group from generators
        %
        % Args:
        %   generators (row cell array of elements of this group): List of generators
        %   order (vpi, optional): Subgroup order
            error('Not implemented');
        end
        
    end
    
end
