classdef SemidirectProductGroup < replab.Group
% Describes an external semidirect product of compact groups
%
% This is an abstract base class.
% Call `.CompactGroup.semidirectProduct` or `.make` to construct an instance.
%
% As semidirect product groups are used as a base for wreath product groups,
% the constructors are duplicated in subclasses as to keep a simple hierarchy of
% constructor calls.
%
% Example:
%   >>> N = replab.S(3);
%   >>> H = replab.PermutationGroup.cyclic(3);
%   >>> A = N.innerAutomorphism([2 3 1]);
%   >>> phi = H.morphismByImages(replab.AutomorphismGroup(N), 'images', {A});
%   >>> sd = H.semidirectProduct(N, @(h, n) phi.imageElement(h).imageElement(n));
%   >>> sd.laws.checkSilent;

    properties (SetAccess = protected)
        H % (`+replab.CompactGroup`): Group acting
        N % (`+replab.CompactGroup`): Group acted upon
        phi % (`+replab.Action`): Action of H on N
    end

    methods (Static) % SemidirectProductGroup creation

        function prd = make(phi)
        % Constructs a semidirect product group from an action
        %
        % Args:
        %   phi (`+replab.Action`): Action of a compact group on another compact group
        %
        % Returns:
        %   `.SemidirectProductGroup`: A specialized instance of `.SemidirectProductGroup`
            isFinite = isa(phi.G, 'replab.FiniteGroup') && isa(phi.P, 'replab.FiniteGroup');
            if isFinite
                prd = replab.prods.SemidirectProductGroup_finite(phi, 'self');
            else
                prd = replab.prods.SemidirectProductGroup_compact(phi);
            end
        end

    end

    methods % Implementations

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
