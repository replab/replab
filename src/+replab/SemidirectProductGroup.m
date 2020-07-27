classdef SemidirectProductGroup < replab.Group
% Describes an external semidirect product of groups
%
% This is an abstract base class. Use `.SemidirectProductGroup.make` or `.CompactGroup.semidirectProduct` to
% construct an instance.

    properties (SetAccess = protected)
        H % (`+replab.Group`): Group acting
        N % (`+replab.Group`): Group acted upon
        phi % (`+replab.Action`): Action of H on N
    end

    methods (Static) % SemidirectProductGroup creation

        function make(phi)
        % Constructs a semidirect product group from an action
        %
        % Args:
        %   phi (`+replab.Action`): Action of a group on another group
        %
        % Returns:
        %   `.SemidirectProductGroup`: A specialized instance of `.SemidirectProductGroup`
            isFinite = all(cellfun(@(g) isa(g, 'replab.FiniteGroup'), factors));
            if isFinite
                prd = replab.prods.SemidirectProductOfFiniteGroups(factors);
            else
                prd = replab.prods.SemidirectProductOfCompactGroups(factors);
            end
        end

    end

    methods

        function self = OfCompactGroups(phi)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            self.phi = phi;
            self.H = H;
            self.N = N;
            self.identity = {H.identity N.identity};
        end

    end

    methods % Implementatoins

        % Domain

        function b = eqv(self, x, y)
            b = self.H.eqv(x{1}, y{1}) && self.N.eqv(x{2}, y{2});
        end

        function g = sample(self)
            g = {self.H.sample self.N.sample};
        end

        % Monoid

        function z = compose(self, x, y)
            xh = x{1};
            xn = x{2};
            yh = y{1};
            yn = y{2};
            % Relation to phi is the conjugation
            % phi_h(n) = h n h^-1
            % we have z = xh xn yh yn = xh yh yh^-1 xn yh yn =
            % = xh yh phi_(yh^-1)(xn) yn
            % and thus
            % zh = xh yh
            % zn = phi_(yh^-1)(xn) yn
            yhinv = self.H.inverse(yh);
            zh = self.H.compose(xh, yh);
            zn = self.N.compose(self.phi.leftAction(yhinv, xn), yn);
            z = {zh zn};
        end

        % Group

        function z = inverse(self, x)
            xh = x{1};
            xn = x{2};
            zh = self.H.inverse(xh);
            zn = self.N.inverse(self.phi.leftAction(xh, xn));
            z = {zh zn};
        end

    end

end
