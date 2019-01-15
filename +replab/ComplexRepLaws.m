classdef ComplexRepLaws < replab.Laws
    properties
        rep;
        G; % group of which rep is a representation
        C; % GL_d(C)
    end
    methods
        function self = ComplexRepLaws(rep)
            self.rep = rep;
            d = self.rep.dimension;
            self.G = rep.group;
            eqvFun = @(X, Y) ~replab.isNonZeroMatrix(X - Y, replab.Settings.doubleEigTol);
            self.C = replab.GroupFun('Mat', eqvFun, @() rand(d, d)+1i*rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
            
        end
        function morphismLaws = laws_asGroupHomomorphism(self)
            morphismLaws = replab.GroupMorphismLaws(@(g) self.rep.image(g), self.G, self.C);
        end
    end
end
