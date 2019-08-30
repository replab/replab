classdef NUSubRep < replab.NURep
% Describes a subrepresentation of a unitary finite representation
    
    properties (SetAccess = protected)
        parent % Parent representation
        F      % Map from the parent to this representation
        G      % Map from this representation to the parent representation
    end
    
    methods
        
        function self = SubRep(parent, F, G)
            d = size(F, 1);
            dParent = size(F, 2);
            assert(size(G, 1) == dParent);
            assert(size(G, 2) == d);
            assert(parent.dimension == dParent);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.F = F;
            self.G = G;
        end

        function s = headerStr(self)
            s = 'Subrepresentation';
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % The projector is expressed in the parent representation
            P = self.G*self.F;
        end
        
        function newSub = collapseParent(self)
        % Returns this subrepresentation as expressed in self.parent.parent
            assert(isa(self.parent, 'replab.NUSubRep'));
            % self W, parent V, parent.parent U
            % self.F: V -> W
            % self.G: W -> V
            % self.parent.F: U -> V
            % self.parent.G: V -> U
            % newF: U -> W
            % newG: W -> U
            newF = self.F * self.parent.F;
            newG = self.parent.G * self.G;
            newSub = self.parent.parent.subRep(newF, newG);
        end
        
        % TODO function sub = in(self, newParent)
        
        % Rep
        
        function rho = image(self, g)
        % Returns G rho(g) F
            rho = self.G * self.parent.image(g) * self.F;
        end
        
    end

end
