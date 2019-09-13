classdef EquivariantDispatch < replab.Dispatch
    
    methods

        function self = EquivariantDispatch
            self.register('ForRepByImages', @(repR, repC) replab.equivariant.ForRepByImages(repR, repC), 10);
            self.register('ForFiniteGroup', @(repR, repC) replab.equivariant.ForFiniteGroup(repR, repC), 5);
            % Default method, works for all compact groups
            self.register('ForCompactGroup', @(repR, repC) replab.equivariant.ForCompactGroup(repR, repC), 0);
        end
        
    end

    methods (Static)
        
        function i = instance
            persistent Instance;
            if isempty(Instance)
                Instance = replab.EquivariantDispatch;
            end
            i = Instance;
        end
        
    end

end
