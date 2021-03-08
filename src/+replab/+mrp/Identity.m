classdef Identity < replab.Isomorphism
% Identity isomorphism

    methods

        function self = Identity(group)
            self.source = group;
            self.target = group;
            if isa(group, 'replab.TorusGroup')
                self.torusMap = eye(group.n);
            elseif isa(group, 'replab.CompactGroup') && group.hasReconstruction
                n = group.maximalTorusDimension;
                self.torusMap = eye(n);
            end
        end

    end

    methods % Implementations

        function y = imageElement(self, x)
            y = x;
        end

        function x = preimageElement(self, y)
            x = y;
        end

    end

    methods (Access = protected)

        function res = computeInverse(self)
            res = self;
        end

    end

end
