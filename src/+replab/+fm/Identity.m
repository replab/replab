classdef Identity < replab.FiniteIsomorphism

    methods

        function self = Identity(group)
            self.source = group;
            self.target = group;
        end

        function x = image(self, x)
        % trivial
        end

        function x = preimage(self, x)
        % trivial
        end

    end

end
