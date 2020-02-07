classdef SubRepNU < replab.Rep
% Describes a subrepresentation of a unitary finite representation

    properties (SetAccess = protected)
        parent % Parent representation
        F      % Map from the parent to this representation
        H      % Map from this representation to the parent representation
    end

    methods

        function self = SubRepNU(parent, F, H)
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

        function [newRep A Ainv] = unitarize(self)
            if isequal(self.parent.isUnitary, true)
                X = self.F * self.F';
                A = chol(X, 'lower');
                U = inv(A) * self.F;
                newRep = self.parent.subRepUnitary(U);
            else
                [newRep A Ainv] = unitarize@replab.Rep(self);
            end
        end

    end

end
