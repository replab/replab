classdef NURep < replab.Str
% Describes a real or complex finite dimension representation of a compact group
%
% EXPERIMENTAL
    
    properties (SetAccess = protected)
        group % Group being representation
        field % 'R' for a representation on a real vector space, 'C' for a representation on a complex vector space
        dimension % Representation dimension
    end
    
    methods % Abstract methods
        
        function rho = image(self, g)
        % Returns the image of a group element (abstract)
        %
        % Args:
        %   g (element of `group`): Element being represented 
        %
        % Returns: 
        %   double matrix: Image of the given element for this representation
            error('Not implemented');
        end
        
    end
    
    methods
        
        function rho = invImage(self, g)
        % Returns the image of the inverse of a group element
        %
        % Args:
        %   g (element of `group`): Element of which the inverse is represented
        %
        % Returns:
        %   double matrix: Image of the inverse of the given element for this representation
            gInv = self.group.inverse(g);
            rho = self.image(gInv);
        end
        
        function b = overR(self)
        % Returns true if this representation is defined over the real field
            b = isequal(self.field, 'R');
        end
        
        function b = overC(self)
        % Returns true if this representation is defined over the complex field
            b = isequal(self.field, 'C');
        end
        
        function rho = sample(self)
            rho = self.image(self.group.sample);
        end
        
        function complexRep = complexification(self)
        % Returns the complexification of a real representation
            assert(self.overR, 'Representation should be real to start with');
            if self.overR
                cRep = self;
            else
                cRep = replab.NURep.lambda(self.group, 'C', self.dimension, @(g) self.image(g));
            end
        end
        
        function sub = subRep(self, F, G)
        % Returns a subrepresentation of this representation
        %
        % Let V be the vector space of this representation, 
        % and W an invariant subspace of V. Let ``F: V -> W``
        % and ``G: W -> V`` be maps such that ``F G = id``,
        % and ``G rho(g) F`` is the subrepresentation.
        % 
        % Args:
        %   F (double matrix): map V -> W
        %   G (double matrix): map W -> V
        % Returns:
        %   NUSubRep: Subrepresentation
            sub = replab.NCSubRep(self, F, G);
        end
        
        %function X1 = actionDualDual(self, g, X)
        %    rho = self.image(g);
        %    rhoi = self.invImage(g);
        %    X1 = rhoi' * X * rho';
        %end
        
        function X1 = innerProductAction(self, g, X)
            rhoi = self.invImage(g);
            X1 = rhoi' * X * rhoi;
        end
                
        function X1 = innerProductProject(self, X)
        % Args:
        %   X (double): Map from V to V*
            T = self.group.decomposition.transversals;
            X1 = X;
            for i = 1:length(T):-1:1
                Ti = T{i};
                S = X1;
                for j = 2:length(Ti)
                    S = S + self.innerProductAction(Ti{j}, X1);
                end
                X1 = S / length(Ti);
            end
        end
        
    end
    
    methods (Static)
        
        function rep = lambda(group, field, dimension, imageFun)
        % Creates a non unitary representation from an image function
            rep = replab.lambda.NURep(group, field, dimension, imageFun);
        end
        
    end
    
end
