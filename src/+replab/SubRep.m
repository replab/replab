classdef SubRep < replab.Rep
% Describes a subrepresentation of a unitary finite representation
    
    properties (SetAccess = protected)
        parent; % Parent representation
                %
                % The basis of this subrepresentation in the parent
                % representation is formed by the rows of the
                % matrix U = diag(D0) * U0, as described below
                %
                % We have image(g) = U * parent.image(g) * U'
                %
        U0;     % Orthogonal but not necessarily normalized basis
                % of the subrepresentation, given in a matrix of
                % dimension dChild x dParent
                %
        D0;     % row vector of correction factors of dimension
                % 1 x dChild
    end
    
    methods
        
        function self = SubRep(parent, U0)
        % Constructs a subrepresentation of the 'parent' representation
        % given by the basis 'U0', which has dimension dThisChild x dParent
            d = size(U0, 1);
            dParent = size(U0, 2);
            assert(parent.dimension == dParent);
            self.group = parent.group;
            self.field = parent.field;
            self.dimension = d;
            self.parent = parent;
            self.isUnitary = true;
            self.U0 = U0;
            D0 = ones(1, d);
            hasCorrection = false;
            for i = 1:d
                nrm = self.U0(i,:) * self.U0(i,:)';
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

        function s = headerStr(self)
            s = 'Subrepresentation';
        end

        function b = hasCorrection(self)
        % Returns true if the basis is not normalized and needs correction
            b = ~isequal(self.D0, []);
        end
        
        function U = U(self)
        % Returns the basis of this subrepresentation in its parent
            if self.hasCorrection
                U = replab.diag_(self.D0) * self.U0;
            else
                U = self.U0;
            end
        end
        
        function sub1 = nice(self)
        % Tries to recover a nice basis for this subrepresentation
            sub1 = replab.rep.niceRep(self);
            if isempty(sub1)
                sub1 = self; % fallback
            end
        end
        
        function P = projector(self)
        % Returns the projector on this subrepresentation
        %
        % The projector is expressed in the parent representation
            P = self.U'*self.U;
        end
        
        function newSub = collapseParent(self)
        % Returns this subrepresentation as expressed in self.parent.parent
            assert(isa(self.parent, 'replab.SubRep'));
            newU0 = self.U0 * self.parent.U;
            newSub = self.parent.parent.subRepUnitary(newU0);
        end
        
        function sub = in(self, newParent)
        % Returns this subrepresentation as expressed in the given 'newParent'
        %
        % newParent must be equal to some self.parent. ... .parent
            if self.parent == newParent % we want handle equality
                sub = self;
            else
                newU0 = self.U0 * self.parent.U;
                newRep = self.parent.parent.subRepUnitary(newU0);
                sub = newRep.in(newParent);
            end
        end
        
        %% Str methods
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            if self.hasCorrection
                for i = 1:size(self.U0, 1)
                    v = replab.str.Normalized(self.U0(i, :), self.D0(i));
                    names{1, end+1} = sprintf('U(%d,:)', i);
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
        
        %% Rep methods
        
        function rho = image(self, g)
            rho = self.U * self.parent.image(g) * self.U';
        end
        
        function rho = inverseImage(self, g)
            rho = self.U * self.parent.inverseImage(g) * self.U';
        end
        
    end

end
