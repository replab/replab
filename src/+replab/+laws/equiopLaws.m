classdef equiopLaws < replab.Laws

    properties (SetAccess = protected)
        op % (`+replab.equiop)
        G % (`+replab.CompactGroup`): Equivariant operator group
        S % (`+replab.Equivariant`): Source equivariant space
        T % (`+replab.Equivariant`): Target equivariant space
    end

    methods

        function self = equiopLaws(op)
            self.op = op;
            self.G = op.group;
            self.S = op.source;
            self.T = op.target;
        end

        function law_equivariant_space_map_S(self, s)
            f = self.op;
            t = double(f(s));
            t1 = self.T.project(t);
            tol = replab.globals.doubleEigTol;
            self.assertApproxEqual(t, t1, tol);
        end

        function law_map_is_equivariant_GS(self, g, s)
            gs = self.S.repR.matrixRowAction(g, self.S.repC.matrixColAction(g, s));
            f = self.op;
            t = double(f(s));
            gt = self.T.repR.matrixRowAction(g, self.T.repC.matrixColAction(g, t));
            gt1 = double(f(gs));
            tol = replab.globals.doubleEigTol;
            self.assertApproxEqual(gt, gt1, tol);
        end

    end

end
