classdef FiniteIdentity < replab.FiniteIsomorphism & replab.fm.Identity
% Identity isomorphism

    methods

        function self = FiniteIdentity(group)
            self@replab.fm.Identity(group);
        end

    end

    methods % Implementations

        function y = imageElement(self, x)
            y = x;
        end

        function x = preimageElement(self, y)
            x = y;
        end

        function T = imageGroup(self, S)
            T = S;
        end

        function S = preimageGroup(self, T)
            S = T;
        end

    end

    methods (Access = protected) % Implementations

        function inv = computeInverse(self)
            inv = self;
        end

    end

end
