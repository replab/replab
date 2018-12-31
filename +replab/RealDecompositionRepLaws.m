classdef RealDecompositionRepLaws < replab.RealRepLaws
%    properties
%        rep;
%        G; % group of which rep is a representation
%        C; % centralizer algebra
%        R; % GL_d(R)
%    end
    methods
        function self = RealDecompositionRepLaws(rep)
            assert(isa(rep, 'replab.RealDecompositionRep'));
            self = self@replab.RealRepLaws(rep);
        end
        function law_parent_orthogonal_subspaces(self)
            tol = replab.Settings.doubleEigTol;
            n = self.rep.nComponents;
            U = cell(1, n);
            for i = 1:n
                c = self.rep.component(i);
                if isempty(c.U)
                    return
                end
                U{i} = c.U;
            end
            for i = 1:n
                for j = i+1:n
                    self.assert(~replab.isNonZeroMatrix(U{i}'*U{j}, tol));
                end
            end
        end
    end
end
