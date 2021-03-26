classdef equiopLaws < replab.Laws

    properties (SetAccess = protected)
        op % (`+replab.equiop)
        S % (`+replab.Equivariant`): Source equivariant space
        T % (`+replab.Equivariant`): Target equivariant space
    end

    methods

        function self = equiopLaws(op)
            self.op = op;
            self.S = op.source;
            self.T = op.target;
        end

        function law_equivariant_space_map_S(self, s)
            f = self.op;
            t = value(f(s));
            t1 = self.T.project(t);
            tol = replab.globals.doubleEigTol;
            self.assertApproxEqual(t, t1, tol);
        end

    end

end
