classdef Action < replab.Obj
% A group action describing the action of a group on a set
%
% Elements of the group ``G`` act upon elements of type ``P``.

    properties (SetAccess = protected)
        G % replab.Group: Group acting
        P % replab.Domain: Set acted upon
    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.ActionLaws(self);
        end

    end

    methods % Action methods

        % Abstract method

        function p1 = leftAction(self, g, p)
        % Computes the left action of a group element on a set element
        %
        % Returns the left action, often denoted by ``p1 = g(p)`` of G over P,
        % which is compatible with group composition as in:
        %
        % `` p2 = g(h(p)) = (g compose h)(p) ``
        %
        % Args:
        %   g (element of `G`): Group element acting
        %   p (element of `P`): Element acted upon
        %
        % Returns:
        %   element of `P`: Result of the left action
            error('Abstract');
        end

        % Default implementations

        function p1 = rightAction(self, p, g)
        % Computes the right action of a group element on a set element
        %
        % Returns the right action, often denoted ``p1 = p^g``,
        % compatible with the group composition as in
        %
        % ``p2 = (p^g)^h = p^(g compose h)``
        %
        % Args:
        %   p (element of `P`): Element acted upon
        %   g (element of `G`): Group element acting
        %
        % Returns:
        %   element of `P`: Result of the right action
            p1 = self.leftAction(self.G.inverse(g), p);
        end

    end

    methods (Static)

        function action = lambda(header, G, P, leftActionFun)
        % Constructs an action from a function handle
        %
        % Args:
        %   G (`+replab.Group`): Group acting
        %   P (`+replab.Domain`: Set acted upon
        %   leftActionFun (function_handle): Handle implementing ``leftAction``
        %
        % Returns:
        %   replab.Action: The action
            action = replab.lambda.Action(header, G, P, leftActionFun);
        end

    end

end
