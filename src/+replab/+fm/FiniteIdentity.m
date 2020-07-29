classdef FiniteIdentity < replab.FiniteIsomorphism & replab.fm.Identity
% Identity isomorphism

    methods

        function self = FiniteIdentity(group)
            self@replab.fm.Identity(group);
        end

    end

    methods % Implementations

        function T = imageGroup(self, S)
            T = S;
        end

        function S = preimageGroup(self, T)
            S = T;
        end

    end

end
