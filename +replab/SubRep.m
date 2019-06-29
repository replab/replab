classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
%
% 
    properties (SetAccess = protected)
        parent; % Parent representation
        U;      % Basis of this subrepresentation in parent
                % U has dimension dThisChild x dParent
        D;      % diagonal matrix of correction factors of dimension
                % dThisChild x dThisChild
                % such that image(g) = D * U * parent.image(g) * U' * D
    end
    
    methods
        
        function self = SubRep(parent, U)
        % U has dimension dThisChild x dParent
            d = size(U, 1);
            dParent = size(U, 2);
            assert(parent.dimension == dParent);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.U = U;
            self.D = speye(d);
            for i = 1:d
                nrm = self.U(i,:) * self.U(i,:)';
                if abs(nrm - 1) > replab.Settings.doubleEigTol
                    self.D(i,i) = 1/sqrt(nrm);
                end
            end
        end
        
        function rho = image(self, g)
            rho = self.D * self.U * self.parent.image(g) * self.U' * self.D;
        end
        
        % TODO optimize action and adjointAction
    end

end
