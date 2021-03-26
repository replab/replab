classdef StandardTorusGroup < replab.TorusGroup
% Describes the standard torus

    methods

        function self = StandardTorusGroup(n)
            self@replab.TorusGroup(zeros(0, n));
        end

    end

end
