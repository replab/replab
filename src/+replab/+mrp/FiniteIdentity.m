classdef FiniteIdentity < replab.FiniteIsomorphism & replab.mrp.Identity
% Identity isomorphism

    methods

        function self = FiniteIdentity(group)
            self@replab.mrp.Identity(group);
        end

    end

    methods % Implementations

        % Morphism

        function y = imageElement(self, x)
            y = x;
        end

        % Isomorphism

        function x = preimageElement(self, y)
            x = y;
        end

        % FiniteMorphism

        function T = imageGroup(self, S)
            T = S;
        end

        function S = preimageGroup(self, T)
            S = T;
        end

        % FiniteIsomorphism

        function l = preservesTypeOrder(self)
            l = true;
        end

    end

    methods (Access = protected) % Implementations

        function inv = computeInverse(self)
            inv = self;
        end

    end

end
