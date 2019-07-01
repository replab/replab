classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
%
% 
    properties (SetAccess = protected)
        parent; % Parent representation
                %
                % The basis of this subrepresentation in the parent
                % representation is formed by the columns of the
                % matrix U = U0 * diag(D0), as described below
                %
                % We have image(g) = U' * parent.image(g) * U
                %
        U0;     % Orthogonal but not necessarily normalized basis
                % of the subrepresentation, given in a matrix of
                % dimension dParent x dChild
                %
        D0;     % row vector of correction factors of dimension
                % 1 x dChild
    end
    
    methods
        
        function self = SubRep(parent, U0)
        % U has dimension dParent x dThisChild
            dParent = size(U0, 1);
            d = size(U0, 2);
            assert(parent.dimension == dParent);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.U0 = U0;
            D0 = ones(1, d);
            hasCorrection = false;
            for i = 1:d
                nrm = self.U0(:,i)' * self.U0(:,i);
                if abs(nrm - 1) > replab.Settings.doubleEigTol
                    hasCorrection = true;
                    D0(i) = 1/sqrt(nrm);
                end
            end
            if hasCorrection
                self.D0 = D0;
            else
                self.D0 = [];
            end
        end
        
        function b = hasCorrection(self)
        % Returns true if the basis is not normalized and needs correction
            b = ~isequal(self.D0, []);
        end
        
        function U = U(self)
            if self.hasCorrection
                U = self.U0 * diag(self.D0);
            else
                U = self.U0;
            end
        end
        
        % Str
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.hasCorrection
                for i = 1:size(self.U0, 2)
                    v = replab.str.Normalized(self.U0(:, i), self.D0(i));
                    names{1, end+1} = sprintf('U(:,%d)', i);
                    values{1, end+1} = v;
                end
            else
                names{1, end+1} = 'U''';
                values{1, end+1} = self.U';
            end
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Rep(self), {'U0' 'D0'} ...
                );
        end
        
        % Rep
        
        function rho = image(self, g)
            rho = self.U' * self.parent.image(g) * self.U;
        end
        
        % TODO optimize action and adjointAction
        
    end

end
