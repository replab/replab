classdef Identity < replab.FiniteIsomorphism
% Identity isomorphism

    methods

        function self = Identity(group)
            self.source = group;
            self.target = group;
        end

        function x = imageElement(self, x)
        % trivial
        end

        function S = imageGroup(self, S)
        % trivial
        end

        function x = preimageElement(self, x)
        % trivial
        end

        function S = preimageGroup(S)
        % trivial
        end

        function res = compose(self, applyFirst)
            res = applyFirst;
        end

        function res = andThen(self, applyLast)
            res = applyLast;
        end

        function res = inverse(self)
            res = self;
        end

    end

end
