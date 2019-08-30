classdef SubRep < replab.nu.Rep
% Describes a subrepresentation of a unitary finite representation
    
    properties (SetAccess = protected)
        parent % Parent representation
        F      % Map from the parent to this representation
        H      % Map from this representation to the parent representation
    end
    
    methods
        
        function self = SubRep(parent, F, H)
            d = size(F, 1);
            dParent = size(F, 2);
            assert(size(H, 1) == dParent);
            assert(size(H, 2) == d);
            assert(parent.dimension == dParent);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.F = F;
            self.H = H;
        end

        function s = headerStr(self)
            s = 'Subrepresentation';
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % The projector is expressed in the parent representation
            P = self.H*self.F;
        end
        
        function newSub = collapseParent(self)
        % Returns this subrepresentation as expressed in self.parent.parent
            assert(isa(self.parent, 'replab.NUSubRep'));
            % self W, parent V, parent.parent U
            % self.F: V -> W
            % self.H: W -> V
            % self.parent.F: U -> V
            % self.parent.H: V -> U
            % newF: U -> W
            % newH: W -> U
            newF = self.F * self.parent.F;
            newH = self.parent.H * self.H;
            newSub = self.parent.parent.subRep(newF, newH);
        end
        
        % TODO function sub = in(self, newParent)
        
        % Rep
        
        function rho = image(self, g)
        % Returns H rho(g) F
            rho = self.H * self.parent.image(g) * self.F;
        end
        
    end

end
